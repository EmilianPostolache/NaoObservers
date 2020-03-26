clear;clc;close;figure;

animated_plot = true;

%% Define walk

foot_distance_x = 0.07;
foot_distance_y = 0.09;
S = 30;                     %% single phase support duration
D = 20;                     %% double phase support duration
N = 100;                    %% called C on paper
Np = 100;                   %% Boh!!
omega = sqrt(9.8/0.33);     %% called ? on paper

balance = false;

if balance == true
    fs_matrix = zeros (26,2);
else
    fs_matrix = [0,-foot_distance_y;%middle of rhe foot point
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
delta = 0.01;       %%Sampling rate used for the MPC horizon
w = 0.05/2;         %%Since in the paper fx, fy (dimension of a rectangular region approximating the footprint) are 0.05, this should be the half of those areas used for dynamic constraints. 

x(1) = 0.0;         %% x coordinate of the CoM
xd(1) = 0.0;        %% derivative of the x coordinate of the CoM
zx(1) = 0.0;        %% x coordinate of the ZMP

y(1) = 0;           %% y coordinate of the CoM
yd(1) = 0;          %% derivative of the y coordinate of the CoM
zy(1) = 0;          %% y coordinate of the ZMP

%% Compute constraints
f1_y = -foot_distance_y;
f2_y = foot_distance_y;

additionalFirstStepDuration = 20;

fs_sequence_x = zeros(S+D+additionalFirstStepDuration,1);       %% Footstep sequence maybe?
fs_sequence_y = zeros(S+D+additionalFirstStepDuration,1);       %% Footstep sequence maybe?

for i = 1:size(fs_matrix,1)-2
    f1_x = fs_matrix(i,1);
    f2_x = fs_matrix(i+1,1);
    
    f1_y = fs_matrix(i,2);
    f2_y = fs_matrix(i+1,2);
    
    fs_sequence_x = [fs_sequence_x; ones(S,1) * f1_x; f1_x + (1:D)'*(f2_x-f1_x)/D];     %%Probably yes since it holds f1 values for S and then linearly interpolate to f2 values in D
    fs_sequence_y = [fs_sequence_y; ones(S,1) * f1_y; f1_y + (1:D)'*(f2_y-f1_y)/D];
end

fs_sequence_x(1) = [];
fs_sequence_y(1) = [];

zx_min = fs_sequence_x - w;         %%Boundary values of the footsteps
zx_max = fs_sequence_x + w;
if balance == true
    zy_min = fs_sequence_y - 5*w;
    zy_max = fs_sequence_y + 5*w;
else
    zy_min = fs_sequence_y - w;
    zy_max = fs_sequence_y + w;
end
zy_min(1:S+D+additionalFirstStepDuration-1) = zy_min(1:S+D+additionalFirstStepDuration-1) - foot_distance_y;
zy_max(1:S+D+additionalFirstStepDuration-1) = zy_max(1:S+D+additionalFirstStepDuration-1) + foot_distance_y;


%% Compute matrices
p = ones(N,1);                      %% 
P = delta*tril(ones(N,N));          %%  No ideas here...  
A = [P;-P];                         %%

%% Compute stability constraint
Aeq = (1-exp(-omega*delta))/omega * exp(-omega*delta*(0:N-1)) - exp(-omega*delta*N) * delta * ones(1,N);


%% CoM constraints
ch = cosh(omega*delta);
sh = sinh(omega*delta);
A_upd = [ch, sh/omega, 1-ch; omega*sh, ch, -omega*sh; 0, 0, 1];
B_upd = [delta-sh/omega; 1-ch; delta];


for i = 1:200  %% ===============>>> PAY ATTENTION HERE!! #iterations
    
    
    P = delta*tril(ones(N,N));
    H = 2*[eye(N)];
    f_x = zeros(1,N);
    f_y = zeros(1,N);
    
    
    %% Disturbances
    
    w_x(i) = 0.09;%%0.08;
    w_y(i) = -0.06;%%0.08;

    
    b_x = [ zx_max(i:i+N-1) - zx(i); - zx_min(i:i+N-1) + zx(i)];
    b_y = [ zy_max(i:i+N-1) - zy(i); - zy_min(i:i+N-1) + zy(i)];
    
    beq_x = x(i) + xd(i)/omega - (1-exp(-omega*N*delta))*zx(i)...
        - omega*delta*exp(-omega*N*delta)*exp(-omega*delta*(0:Np-1))*fs_sequence_x(i+N:i+N+Np-1)...
        - exp(-omega*(N+Np)*delta)*fs_sequence_x(i+(N+Np));
    beq_y = y(i) + yd(i)/omega - (1-exp(-omega*N*delta))*zy(i)...
        - omega*delta*exp(-omega*N*delta)*exp(-omega*delta*(0:Np-1))*fs_sequence_y(i+N:i+N+Np-1)...
        - exp(-omega*(N+Np)*delta)*fs_sequence_y(i+(N+Np));
    
    zd_x = quadprog(H,f_x,A,b_x,Aeq,beq_x);
    zd_y = quadprog(H,f_y,A,b_y,Aeq,beq_y);
    
    
    z_pred_x = P*zd_x + zx(i);
    z_pred_y = P*zd_y + zy(i);
    
    
    
    x_updated = A_upd*[x(i); xd(i); zx(i)] + B_upd*(zd_x(1))+[0;delta;0]*w_x(i);
    y_updated = A_upd*[y(i); yd(i); zy(i)] + B_upd*(zd_y(1))+[0;delta;0]*w_y(i);
    
    
    x(i+1) = x_updated(1);
    xd(i+1) = x_updated(2);
    zx(i+1) = x_updated(3);
    
    y(i+1) = y_updated(1);
    yd(i+1) = y_updated(2);
    zy(i+1) = y_updated(3);
    
    %% Divergent component of motion
    
    x_u(i) = x(1,i)+xd(1,i)/omega;
    y_u(i) = y(1,i)+yd(1,i)/omega;
    
    %% XY Plot
    
    if (animated_plot)
        
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
    
    
    disp('iteration number ')
    disp(i)
    
end


%% XY Plot

if (~animated_plot)
    
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