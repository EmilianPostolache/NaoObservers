clear;clc;close;figure;

animated_plot = true;

%% Define walk

foot_distance_x = 0.07;     %% distance between centroids of two consecutive footsteps (along x axis) [m]
foot_distance_y = 0.09;     %% distance between centroids of two consecutive footsteps (along y axis) [m]
S = 30;                     %% single phase support duration
D = 20;                     %% double phase support duration
N = 100;                    %% called C on paper (size control horizon)
Np = 100;                   %% Boh!!
omega = sqrt(9.8/0.33);     %% called eta on paper (0.33 is the constant height of the CoM)[m]

balance = false;

if balance == true
    fs_matrix = zeros (26,2);              % 26 is the number of footsteps planned
else
    fs_matrix = [0,-foot_distance_y;       %middle of the foot point (centroid step)
        foot_distance_x,foot_distance_y;
        2*foot_distance_x,-foot_distance_y;
        3*foot_distance_x,foot_distance_y;
        4*foot_distance_x,-foot_distance_y;
        5*foot_distance_x,foot_distance_y;
        6*foot_distance_x,-foot_distance_y;
        7*foot_distance_x,foot_distance_y;
        8*foot_distance_x,-foot_distance_y;
        9*foot_distance_x,+foot_distance_y;
        10*foot_distance_x,-foot_distance_y;
        11*foot_distance_x,+foot_distance_y;
        12*foot_distance_x,-foot_distance_y
        13*foot_distance_x,+foot_distance_y;
        14*foot_distance_x,-foot_distance_y;
        15*foot_distance_x,+foot_distance_y
        16*foot_distance_x,-foot_distance_y;
        17*foot_distance_x,+foot_distance_y;
        18*foot_distance_x,-foot_distance_y;
        19*foot_distance_x,+foot_distance_y;
        20*foot_distance_x,-foot_distance_y;
        21*foot_distance_x,+foot_distance_y;
        22*foot_distance_x,-foot_distance_y;
        23*foot_distance_x,+foot_distance_y;
        24*foot_distance_x,-foot_distance_y;
        25*foot_distance_x,+foot_distance_y;
        26*foot_distance_x,-foot_distance_y];
end

%% General parameters
delta = 0.01;       % duration sampling intervals used for the MPC horizon [s]
w = 0.05/2;         %%Since in the paper fx, fy (dimension of a rectangular region approximating the footprint) are 0.05, this should be the half of those areas used for dynamic constraints. 
                    % margin footsteps
x(1) = 0.0;         % initial value of CoM position (x axis)
xd(1) = 0.0;        % initial value of CoM velocity (x axis)
zx(1) = 0.0;        % initial value of ZMP position (x axis)

y(1) = 0;           % initial value of CoM position (y axis)
yd(1) = 0;          % initial value of CoM velocity (y axis)
zy(1) = 0;          % initial value of ZMP position (y axis)

%% Compute constraints
f1_y = -foot_distance_y;
f2_y = foot_distance_y;

additionalFirstStepDuration = 20;  % gives extra time to stabilize first step cycle (SS+DS)

fs_sequence_x = zeros(S+D+additionalFirstStepDuration,1);       % temporal behaviour of position ZMP (along x asis)
fs_sequence_y = zeros(S+D+additionalFirstStepDuration,1);       % temporal behaviour of position ZMP (along y asis)

for i = 1:size(fs_matrix,1)-2 % for 24 footsteps (2 less than original 26), maybe boundary conditions?
    f1_x = fs_matrix(i,1);   % get current value of ZMP pos x
    f2_x = fs_matrix(i+1,1); % get next value of ZMP pos x

    f1_y = fs_matrix(i,2);   % same for y (actually this value oscillates)
    f2_y = fs_matrix(i+1,2);
    % behaviour ZMP:
    fs_sequence_x = [fs_sequence_x; ones(S,1) * f1_x; f1_x + (1:D)'*(f2_x-f1_x)/D];     % column vector that increases size after each iteration by 50 (initially 70): during SS (duration:20) value constant to current one (f1_x), then linearly interpolation (lasts duration DS) to next value (f2_x). Results in a stair-like graph (steps are slopes).
    fs_sequence_y = [fs_sequence_y; ones(S,1) * f1_y; f1_y + (1:D)'*(f2_y-f1_y)/D];     % same for y but here classical oscillating behaviour 
end

fs_sequence_x(1) = [];
fs_sequence_y(1) = [];

zx_min = fs_sequence_x - w;         % Subtract margin footstep for minimum value for ZMP (x axis)
zx_max = fs_sequence_x + w;         % Add margin footstep for maximum value for ZMP (x axis)
if balance == true                  % for balance=true larger margin for ZMP (fs_matrix = zeros!)
    zy_min = fs_sequence_y - 5*w;
    zy_max = fs_sequence_y + 5*w;
else
    zy_min = fs_sequence_y - w;     % otherwise normal values
    zy_max = fs_sequence_y + w;
end
zy_min(1:S+D+additionalFirstStepDuration-1) = zy_min(1:S+D+additionalFirstStepDuration-1) - foot_distance_y;  % fix first step ZMP boundary conditions (aligned footsteps)
zy_max(1:S+D+additionalFirstStepDuration-1) = zy_max(1:S+D+additionalFirstStepDuration-1) + foot_distance_y;  % same


%% Compute matrices                 % same meaning of page 3 of paper "Intrinsically Stable MPC for Humanoid Gait Generation"
p = ones(N,1);                      % column vector size N
P = delta*tril(ones(N,N));          % extracts lower triangular matrix (obtains a matrix where main diagonal and lower part =1 while upper part =0)
A = [P;-P];                         % get 2 lower triangular matrices stacked (200x100) for ZMP constraints (see notes)

%% Compute stability constraint
Aeq = (1-exp(-omega*delta))/omega * exp(-omega*delta*(0:N-1)) - exp(-omega*delta*N) * delta * ones(1,N);  % parameter for optimization function quadprog (Aeq*x=beq) as constraint


%% CoM constraints
ch = cosh(omega*delta);
sh = sinh(omega*delta);
A_upd = [ch, sh/omega, 1-ch; omega*sh, ch, -omega*sh; 0, 0, 1];  % rotation matrix?
B_upd = [delta-sh/omega; 1-ch; delta];


for i = 1:200  %% ===============>>> PAY ATTENTION HERE!! #iterations

    P = delta*tril(ones(N,N));  % same as row 99
    H = 2*[eye(N)];  % identity matrix with 2 on main diagonal (because quadprog minimize 0.5*H not only H)
    f_x = zeros(1,N);  % row vector, no need of adding elements to H in quadprog
    f_y = zeros(1,N);  % same


    %% Disturbances

    w_x(i) = 0.09;%%0.08;  % disturbances along x axis at instant i [m/s^2]
    w_y(i) = -0.06;%%0.08; % disturbances along y axis at instant i [m/s^2]


    b_x = [ zx_max(i:i+N-1) - zx(i); - zx_min(i:i+N-1) + zx(i)];  % why not zx_max(i+1:i+N)? b part of ZMP contraints with ZMP values extracted from footstep sequence +- safety margin(see notes)
    b_y = [ zy_max(i:i+N-1) - zy(i); - zy_min(i:i+N-1) + zy(i)];  % same for y
0.
    beq_x = x(i) + xd(i)/omega - (1-exp(-omega*N*delta))*zx(i)...
        - omega*delta*exp(-omega*N*delta)*exp(-omega*delta*(0:Np-1))*fs_sequence_x(i+N:i+N+Np-1)...
        - exp(-omega*(N+Np)*delta)*fs_sequence_x(i+(N+Np));  % for quadprog see constraints
    beq_y = y(i) + yd(i)/omega - (1-exp(-omega*N*delta))*zy(i)...
        - omega*delta*exp(-omega*N*delta)*exp(-omega*delta*(0:Np-1))*fs_sequence_y(i+N:i+N+Np-1)...
        - exp(-omega*(N+Np)*delta)*fs_sequence_y(i+(N+Np));
    % decoupled solving
    zd_x = quadprog(H,f_x,A,b_x,Aeq,beq_x);  % solve current quadratic programming problem for x
    zd_y = quadprog(H,f_y,A,b_y,Aeq,beq_y);  % same for y


    z_pred_x = P*zd_x + zx(i);  % formula 7 paper "Intrinsically Stable MPC for Humanoid Gait Generation"
    z_pred_y = P*zd_y + zy(i);  % same (both not used)



    x_updated = A_upd*[x(i); xd(i); zx(i)] + B_upd*(zd_x(1))+[0;delta;0]*w_x(i); % why disturbance only velocity CoM? Apply first value of optimal input sequence obtained from MPC
    y_updated = A_upd*[y(i); yd(i); zy(i)] + B_upd*(zd_y(1))+[0;delta;0]*w_y(i); % update rule from LIP dynamic integration?


    x(i+1) = x_updated(1);  % update components state vector
    xd(i+1) = x_updated(2);
    zx(i+1) = x_updated(3);

    y(i+1) = y_updated(1);
    yd(i+1) = y_updated(2);
    zy(i+1) = y_updated(3);

    %% Divergent component of motion

    x_u(i) = x(1,i)+xd(1,i)/omega;  % unstable subsystem (always decoupled)
    y_u(i) = y(1,i)+yd(1,i)/omega;  % same

    %% XY Plot

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

        h2 = plot(x,y,'r','lineWidth',2);   % plot behaviour CoM (red)
        h3 = plot(zx,zy,'b','lineWidth',2); % plot behaviour ZMP (blue)

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


    disp('iteration number ')
    disp(i)

end


%% XY Plot

if (~animated_plot)  % same

    figure(1)
    clf
    hold on
    grid on
    rect_x = [w,w,-w,-w,w];
    rect_y = [foot_distance_y+w,-foot_distance_y-w,-foot_distance_y-w,foot_distance_y+w,foot_distance_y+w];
    plot(rect_x,rect_y,'m','lineWidth',2,'HandleVisibility','off');

    rect_x = [w,w,-w,-w,w];
    rect_y = [w,-w,-w,w,w];

    nPlottedFootsteps = 14;

    for j = 1:nPlottedFootsteps
        rect_x = [w,w,-w,-w,w];
        rect_y = [w,-w,-w,w,w];
        h1 = plot(fs_matrix(j,1)+rect_x,fs_matrix(j,2)+rect_y,'m','lineWidth',2,'HandleVisibility','off');
    end

    h2 = plot(x,y,'r','lineWidth',2);
    h3 = plot(zx,zy,'b','lineWidth',2);

    legend('CoM', 'ZMP')
    axis equal
    if balance == true
        axis([-0.1 0.1 -0.3 0.3])
    else
        axis([-0.2 1 -0.2 0.2])
    end
    xlabel('x [m]')
    ylabel('y [m]')
end
