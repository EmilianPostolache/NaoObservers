clear;clc;close;figure;

animated_plot = true;

%% Define walk

foot_distance_x = 0.07;        % distance between centroids of two consecutive footsteps (along x axis) [m]
foot_distance_y = 0.09;        % distance between centroids of two consecutive footsteps (along y axis) [m]
n_footsteps     = 26;          % number of footsteps
balance         = false;       % false if walking
g               = 9.8;         % gravity constant
hc              = 0.33;        % (constant height of the CoM)[m]
eta             = sqrt(g/hc);  %
S               = 30;          % single phase support duration
D               = 20;          % double phase support duration
NC              = 100;         % called C on paper (size control horizon)
NP              = 100;         % called P on paper (prediction horizon)
iters           = 200;

if balance == true
    fs_matrix = zeros(n_footsteps, 2);
else
    fs_matrix = zeros(n_footsteps, 2);
    for k=0:n_footsteps-1
        fs_matrix(k+1, 1) = k * foot_distance_x;
        fs_matrix(k+1, 2) = (-1)^(k+1) * foot_distance_y;
    end
end

%% General parameters
delta     = 0.01;              % duration sampling intervals used for the MPC horizon [s]
w         = 0.05/2;            % since in the paper fx, fy (dimension of a rectangular region approximating the footprint) are 0.05, this should be the half of those areas used for dynamic constraints.

xc        = zeros(1, iters);   % CoM (x axis)
xc_dot    = zeros(1, iters);   % CoM velocity (x axis)
xz        = zeros(1, iters);   % ZMP position (x axis)

yc        = zeros(1, iters);   % CoM (y axis)
yc_dot    = zeros(1, iters);   % CoM velocity (y axis)
yz        = zeros(1, iters);   % ZMP position (y axis)

xc(1)     = 0.0;     % initial value of CoM position (x axis)
xz(1)     = 0.0;     % initial value of ZMP position (x axis)
xc_dot(1) = 0.0;     % initial value of CoM velocity (x axis)

yc(1)     = 0;       % initial value of CoM position (y axis)
yc_dot(1) = 0;       % initial value of CoM velocity (y axis)
yz(1)     = 0;       % initial value of ZMP position (y axis)

%% Compute constraints
add_steps = 20;  % gives extra time to stabilize first step cycle (SS+DS)

fs_sequence_x = zeros(add_steps + length(fs_matrix) * (S + D), 1); % temporal behaviour of position ZMP (along x asis)
fs_sequence_y = zeros(add_steps + length(fs_matrix) * (S + D), 1); % temporal behaviour of position ZMP (along y asis)

for k = 1:length(fs_matrix)-1   % for 24 footsteps (2 less than original 26), maybe boundary conditions?
    f1_x = fs_matrix(k,1);      % get current value of ZMP pos x
    f2_x = fs_matrix(k+1,1);    % get next value of ZMP pos x

    f1_y = fs_matrix(k,2);      % same for y (actually this value oscillates)
    f2_y = fs_matrix(k+1,2);
    
    % behaviour ZMP:
    fs_sequence_x(add_steps + k*(S+D) + 1: ...
                  add_steps + k*(S+D) + S, 1) = f1_x;
    fs_sequence_x(add_steps + k*(S+D) + S + 1: ...
                  add_steps + k*(S+D) + S + D, 1) = f1_x + (1:D)' * (f2_x-f1_x)/D;   % column vector that increases size after each iteration by 50 (initially 70): during SS (duration:20) value constant to current one (f1_x), then linearly interpolation (lasts duration DS) to next value (f2_x). Results in a stair-like graph (steps are slopes).
              
    fs_sequence_y(add_steps + k*(S+D) + 1: ...
                  add_steps + k*(S+D) + S, 1) = f1_y;
    fs_sequence_y(add_steps + k*(S+D) + S + 1: ...
                  add_steps + k*(S+D) + S + D, 1) = f1_y + (1:D)' * (f2_y-f1_y)/D;   % same for y but here classical oscillating behaviour         
end

xz_min = fs_sequence_x - w;         % Subtract margin footstep for minimum value for ZMP (x axis)
xz_max = fs_sequence_x + w;         % Add margin footstep for maximum value for ZMP (x axis)
if balance == true                  % for balance=true larger margin for ZMP (fs_matrix = zeros!)
    yz_min = fs_sequence_y - 5*w;
    yz_max = fs_sequence_y + 5*w;
else
    yz_min = fs_sequence_y - w;     % otherwise normal values
    yz_max = fs_sequence_y + w;
end
yz_min(1:S+D+add_steps) = yz_min(1:S+D+add_steps) - foot_distance_y;  % fix first step ZMP boundary conditions (aligned footsteps)
yz_max(1:S+D+add_steps) = yz_max(1:S+D+add_steps) + foot_distance_y;  % same


%% Compute matrices                 
p = ones(NC,1);                     % column vector size NC; same meaning of page 3 of paper "Intrinsically Stable MPC for Humanoid Gait Generation"
P = delta*tril(ones(NC,NC));        % extracts lower triangular matrix (obtains a matrix where main diagonal and lower part =1 while upper part =0)
A = [P;-P];                         % get 2 lower triangular matrices stacked (200x100) for ZMP constraints (see notes)
H = 2*eye(NC);                      % identity matrix with 2 on main diagonal (because quadprog minimize 0.5*H not only H)
f = zeros(1,NC);                    % row vector, no need of adding elements to H in quadprog

%% Compute stability constraint
lambda = exp(-eta*delta);

Aeq = (1-lambda)/eta * exp(-eta*delta*(0:NC-1)) - exp(-eta*delta*NC) * delta * ones(1,NC);  % parameter for optimization function quadprog (Aeq*x=beq) as constraint

%% CoM constraints
ch = cosh(eta*delta);
sh = sinh(eta*delta);
A_upd = [ch, sh/eta, 1-ch; eta*sh, ch, -eta*sh; 0, 0, 1];  % solution to LIP equation (see pag.130 of Introduction to Humanoid Robotics)
B_upd = [delta-sh/eta; 1-ch; delta];                       

%% Disturbances
d_x = 0.09;   % disturbances along x axis at instant k [m/s^2] (stable: 0.08)
d_y = -0.06;  % disturbances along y axis at instant k [m/s^2] (stable: 0.08)

%% Main Loop
for k = 1:iters  % PAY ATTENTION HERE!! #iterations

    b_x = [xz_max(k+1:k+NC) - xz(k); - xz_min(k:k+NC-1) + xz(k)];  % why not xz_max(k+1:k+NC)? b part of ZMP contraints with ZMP values extracted from footstep sequence +- safety margin(see notes)
    b_y = [yz_max(k+1:k+NC) - yz(k); - yz_min(k:k+NC-1) + yz(k)];  % same for y
    
    beq_x = xc(k) + xc_dot(k)/eta - (1-exp(-eta*NC*delta))*xz(k)...
        - eta*delta*exp(-eta*NC*delta)*exp(-eta*delta*(0:NP-1))*fs_sequence_x(k+NC:k+NC+NP-1)...
        - exp(-eta*(NC+NP)*delta)*fs_sequence_x(k+(NC+NP));  % for quadprog see constraints
    beq_y = yc(k) + yc_dot(k)/eta - (1-exp(-eta*NC*delta))*yz(k)...
        - eta*delta*exp(-eta*NC*delta)*exp(-eta*delta*(0:NP-1))*fs_sequence_y(k+NC:k+NC+NP-1)...
        - exp(-eta*(NC+NP)*delta)*fs_sequence_y(k+(NC+NP));

    xz_dot = quadprog(H,f,A,b_x,Aeq,beq_x);  % solve current quadratic programming problem for x
    yz_dot = quadprog(H,f,A,b_y,Aeq,beq_y);  % same for y

    x_updated = A_upd*[xc(k); xc_dot(k); xz(k)] + B_upd*(xz_dot(1))+[0;delta;0]*d_x; % why disturbance only velocity CoM? Apply first value of optimal input sequence obtained from MPC
    y_updated = A_upd*[yc(k); yc_dot(k); yz(k)] + B_upd*(yz_dot(1))+[0;delta;0]*d_y; % update rule from LIP dynamic integration

    xc(k+1) = x_updated(1);  % update components state vector
    xc_dot(k+1) = x_updated(2);
    xz(k+1) = x_updated(3);

    yc(k+1) = y_updated(1);
    yc_dot(k+1) = y_updated(2);
    yz(k+1) = y_updated(3);

    % XY Plot
    if (animated_plot)

        figure(1)
        clf
        hold on    % show axes
        grid on    % show grid
        rect_x = [w,w,-w,-w,w];
        rect_y = [foot_distance_y+w,-foot_distance_y-w,-foot_distance_y-w,foot_distance_y+w,foot_distance_y+w];
        plot(rect_x,rect_y,'m','lineWidth',2,'HandleVisibility','off'); % plot robot initial position (aligned feet)

        rect_x = [w,w,-w,-w,w];  % measure foot
        rect_y = [w,-w,-w,w,w];

        nPlottedFootsteps = 14;  % how many footsteps to plot

        for j = 1:nPlottedFootsteps
            rect_x = [w,w,-w,-w,w];  % same as lines 176-177
            rect_y = [w,-w,-w,w,w];
            h1 = plot(fs_matrix(j,1)+rect_x,fs_matrix(j,2)+rect_y,'m','lineWidth',2,'HandleVisibility','off'); % plot footsteps centered on sequence
        end

        h2 = plot(xc,yc,'r','lineWidth',2);   % plot behaviour CoM (red)
        h3 = plot(xz,yz,'b','lineWidth',2); % plot behaviour ZMP (blue)

        legend('CoM', 'ZMP')
        axis equal
        if balance == true
            axis([-0.1 0.1 -0.3 0.3])  % scale axes
        else
            axis([-0.2 1 -0.2 0.2])    % scale axes
        end
        xlabel('x [m]')
        ylabel('y [m]')

    end
    disp('Iteration number: ')
    disp(k)
end