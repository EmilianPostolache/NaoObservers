%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% F.M. Smaldone 
%
% EXPLORING BOUNDEDNESS CONDITION IN A MPC FRAMEWORK
%
% Problem: ZMP contraints need to be satisfied -->
% --> Not all the disturbances can be bounded.
%
% Question: is it only a matter of feasibility recovery 
% or also a disturbance compensation should be performed?
%
% W.l.o.g. the problem is only considered along x-axis.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc;close;figure;

%% Define walk

foot_distance_x = 0.05;
foot_distance_y = 0.18;
S = 20;
D = 10;
N = 60;
Np = 105;
omega = sqrt(9.8/0.33);
m = 4.5; % Robot Mass

%% Restriction trial
max_dist = 0.4; %max_dist = 0.3 for (true) disturbance
delta = 0.01;
max_vel = max_dist/delta;

r = 2.6210e-04; % exact term
r = 3.0210e-04; % more robust term
r = 5.0210e-04; % more robust term
r_computed = 3*0.4*delta*exp(omega*delta)/omega^2 + 0.4*delta/omega^2;
r_computed = delta*(3*max_dist*delta*exp(omega*delta)/omega^2 + max_dist*(exp(omega*delta)-1)/omega^2)/(1-exp(-omega*0.39));
r_computed = delta*(delta*max_vel/omega^2 + max_dist*(exp(omega*delta)-1)/omega^2)/(1-exp(-omega*N*delta));
r_computed = delta*(4*max_dist*(exp(omega*delta)+1)/omega^2 + max_dist*(exp(omega*delta)-1)/omega^2)/(1-exp(-omega*N*delta));
r_computed = delta*(4*max_dist*(exp(omega*delta)+1)/omega^2 + max_dist*(exp(omega*delta)-1)/omega^2)/(1-exp(-omega*19*delta));

alpha = 2;
r_computed = delta*(2*alpha*max_dist*(exp(omega*delta)+1)/omega^2 + max_dist*(exp(omega*delta)-1)/omega^2)/(1-exp(-omega*13*delta)); % for alpha = 2
r = r_computed;
restriction_x = r*(0:N-1)';
restriction_x(13:50) = restriction_x(12); % for alpha = 2

% alpha = 2.5;
% r_computed = delta*(2*alpha*max_dist*(exp(omega*delta)+1)/omega^2 + max_dist*(exp(omega*delta)-1)/omega^2)/(1-exp(-omega*4*delta)); % for alpha = 2
% r = r_computed;
% restriction_x = r*(0:N-1)';
% restriction_x(5:50) = restriction_x(4); % for alpha = 2


pulse = 5;
ampl = 0.25;

restriction_x = 0*[restriction_x;restriction_x];

%% Choose whether to perform balancing of gait

balance = false;

if balance == true
    fs_matrix = zeros(90,2);
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
delta = 0.01;
w = 0.08/2;
w = 0.025;

x(1) = 0.0;
xd(1) = 0.0;
zx(1) = 0.0;

y(1) = 0;
yd(1) = 0;
zy(1) = 0;

x_u(1) = 0;

%% Compute constraints
f1_y = -foot_distance_y;
f2_y = foot_distance_y;

additionalFirstStepDuration = 0*50;

fs_sequence_x = zeros(S+D+additionalFirstStepDuration,1);
fs_sequence_y = zeros(S+D+additionalFirstStepDuration,1);

for i = 1:26
    f1_x = fs_matrix(i,1);
    f2_x = fs_matrix(i+1,1);
    
    f1_y = fs_matrix(i,2);
    f2_y = fs_matrix(i+1,2);
    
    fs_sequence_x = [fs_sequence_x; ones(S,1) * f1_x; f1_x + (1:D)'*(f2_x-f1_x)/D];
    fs_sequence_y = [fs_sequence_y; ones(S,1) * f1_y; f1_y + (1:D)'*(f2_y-f1_y)/D];
end

fs_sequence_x(1) = [];
fs_sequence_y(1) = [];

zx_min = fs_sequence_x - w;
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
p = ones(N,1);
P = delta*tril(ones(N,N));
A = [P;-P];

%% Compute stability constraint
Aeq = (1-exp(-omega*delta))/omega * exp(-omega*delta*(0:N-1)) - exp(-omega*delta*N) * delta * ones(1,N);

%% Initialize Disturbance Observer
d_z_x(1)= 15;
% x_obs = [x(1);xd(1);zx(1)];
y_obs = [y(1);yd(1);zy(1)];

des_x(1)=0;
des_y(1)=0;
global dx  dy;
dx = 0;
dy = 0;
des_x_mean(1)=0;
des_y_mean(1)=0;

%% CoM constraints
    ch = cosh(omega*delta);
    sh = sinh(omega*delta);
    global A_upd B_upd;
    A_upd = [ch, sh/omega, 1-ch; omega*sh, ch, -omega*sh; 0, 0, 1];
    B_upd = [delta-sh/omega; 1-ch; delta];
%% Matrix that may come in handy
    num_states=3;
    C = [0 0 1];
    
    S_bar = zeros(N, N);
    for j = 1:N
        for k = 1:j
            S_bar(j,k) = C*A_upd^(j-k)*B_upd;
        end
    end
    
    T_bar = zeros(N, 3);
    for j = 1:N
        T_bar(j, :) = C*A_upd^j;
    end
    
    S_bar_state = zeros(num_states*N, N);
    for j = 1:N
        for k = 1:j
            S_bar_state(num_states*(j-1)+1:j*num_states,k) = A_upd^(j-k)*B_upd;
        end
    end
    
    T_bar_state = zeros( num_states*N, num_states);
    for j = 1:N
        T_bar_state(num_states*(j-1)+1:j*num_states,:) = A_upd^j;
    end
    
  T=[-T_bar;T_bar]; %goes with plus! 
  
  A_d=[0,1;omega^2,0];B_d=[0;1];
  sys=ss(A_d,B_d,[1 0],[]);
  sys=c2d(sys,delta);
  A_d=sys.A;
  B_d=sys.B;
  
  X(:,1)=[0;0]; Y(:,1)=[0;0];
  
  for i=1:size(fs_sequence_x,1)
  X(:,i+1)=A_d*X(:,i)+B_d*fs_sequence_x(i,1);
  Y(:,i+1)=A_d*Y(:,i)+B_d*fs_sequence_y(i,1);
  end
 
  D=[0;1;0];
  D_bar = zeros(N*num_states, N);
  
  T_bar_state_D = zeros( N, num_states);
  for j = 1:N
      T_bar_state_D(j,:) = [0 1 0]*A_upd^j;
  end
  D_ZMP = T_bar_state_D*D*(1/omega^2);
  
  
    global C;
    C = [0 1 0];
    D_bar = zeros(N*num_states, N);
  
   S_bar_state = zeros(N, N);
   S_bar_state_D = zeros(N, N);
    for j = 1:N
        for k = 1:j
            S_bar_state(j,k) = C*A_upd^(j-k)*B_upd;
            S_bar_state_D(j,k) = C*A_upd^(j-k)*B_upd;
        end
    end
    
 %% Matrices for MPC
 
      P = delta*tril(ones(N,N));
      H = 2*[eye(N)];
      f_x = zeros(1,N);
      f_y = zeros(1,N);
      
   

      %% Solve: MPC and cloed loop evolution
      Prova = 0;  
         %% Useful vectors
      for i = 0:N-1
          Deltas(i+1) = delta*exp(-omega*delta*i);
      end
      if Prova == 1 
      for i = 0:4*N-1
          exps(i+1) = exp(-omega*delta*i);
      end
      elseif Prova == 4
           for i = 0:N-1
          exps(i+1) = exp(-omega*delta*i);
           end
      elseif Prova == 5  % periodic tail 
           for i = 0:N-1
          exps(i+1) = exp(-omega*delta*i);
           end
      end

      for i = 0:N-1
          exps(i+1) = exp(-omega*delta*i);
      end
      
%% initialize observer 

x_obs(:,1)=[0;0;0;0;0];
A_obs = [0 1 0 0 0; omega^2 0 -omega^2 1 0; 0 0 0 0 0; 0 0 0 0 1; 0 0 0 0 0 ];
B_obs = [0;0;1;0;0];


 p = [-10 -20 -22 -15 -16];
 C_obs = [1 0 0 0 0;0 0 1 0 0];
%   p = [-1 -2 -3 -4 -5];
 L = place((A_obs)',[1 0 0 0 0;0 0 1 0 0]',p);
sys_for_lqr = ss((A_obs)',[1 0 0 0 0;0 0 1 0 0]', [0 0 0 1 0; 1 0 0 0 0], []);
% L = lqr(sys_for_lqr, diag([0.01, 0.0001, 1/0.0250, 1 , 10*(1/0.1814)]), [0.00001 0; 0 0.01],zeros(5,2));
%  L = lqr(sys_for_lqr, diag([0.0001, 0.1, 1/0.0250, 10 , 1000*(1/0.1814)]), [0.00001 0; 0 0.01],zeros(5,2));

 disturbance(1:400) = 0.2;

 
      for i = 1:400
        
          
          
          
      %% Disturbance
      if Prova == 0 || Prova == 3
          A_disturbancex = (delta*0.25);
          A_disturbancex = disturbance(i);
          DDXX = 0*x_obs(4,i);
          F(i) = (A_disturbancex);
          knownDisturbance = F(i);
          
          linear_component = 0*exps*ones(size(exps,2),1)*x_obs(5,i)/omega^3;

          
      elseif Prova == 1
          A_disturbancex = (delta*0.3)*sin(2*pi*i*delta/7);
          F(i) = (A_disturbancex);
          knownDisturbance = F(i);
          for j = i:i+4*N-1
              D_x(j-i+1)= (delta*0.3)*sin(2*pi*j*delta/7);
          end
          integral = exps*D_x';
      elseif Prova == 4
          A_disturbancex = (delta*0.3)*sin(2*pi*i*delta/7);
          F(i) = (A_disturbancex);
          knownDisturbance = F(i);
          for j = i:i+N-1
              D_x(j-i+1)= (delta*0.3)*sin(2*pi*j*delta/7);
          end
          integral = exps*D_x';
      elseif Prova == 5
          A_disturbancex = (delta*0.15)+(delta*0.02)*sin(2*pi*i*delta/7);
%            A_disturbancex = (delta*0.3)*sin(2*pi*i*delta/7);
          F(i) = (A_disturbancex);
          knownDisturbance = F(i);
          for j = i:i+N-1
              D_x(j-i+1)= (delta*0.2)*sin(2*pi*j*delta/7);
          end
          integral = exps*D_x';
          integral2 = 0;
%           NN = int16(N);
%           NNP = int16(Np);
          n_period = fix(Np/N);
          m = mod(Np,N);
          for j = 0:n_period-1
              exps2 = []; D_x2 = [];
          for k = i+(1+j)*N:i+(2+j)*N
              exps2 = [exps2,exp(-omega*delta*(k))];
              D_x2 = [D_x2,(delta*0.2)*sin(2*pi*(k-(1+j)*N)*delta/7)];
          end
          integral2 = integral2 + exps2*D_x2';
          end 
          exps3 = []; D_x3 = [];
          for j = i+n_period*N:i+n_period*N+m
              exps3 = [exps3,exp(-omega*delta*(k))];
              D_x3 = [D_x3,(delta*0.2)*sin(2*pi*(j-n_period*N)*delta/7)];
          end
          integral3 = exps3*D_x3';
          integral = integral + integral2 + integral3 + delta*0.15*(1/(1-exp(-omega*delta)));
      end
%       A_disturbancex = +(delta*0.35);
%       F(i) = (A_disturbancex+0.1*A_disturbancex*sin(2*pi*i/100)+0.1*A_disturbancex*sin(2*pi*i/120));
%       knownDisturbance = F(i) + normrnd(0,0.1*delta*0.35);
      % very nice case with A_disturbancex = +(delta*0.35) or 0.3 and
      % knownDisturbance = F(i) !
      
      %% Further compensation

      
      %% Constraints
      b_x = [ zx_max(i:i+N-1) - zx(i); - zx_min(i:i+N-1) + zx(i)] - restriction_x ; %- restriction_x
      
      if Prova == 0
      beq_x = x(i) + xd(i)/omega - (1-exp(-omega*N*delta))*zx(i)+(1/omega^2)*((DDXX))+0*(1/omega)*(delta*(DDXX))*(1/(1-exp(-omega*delta)))+linear_component/omega^2 ...
          - omega*delta*exp(-omega*N*delta)*exp(-omega*delta*(0:Np-1))*fs_sequence_x(i+N:i+N+Np-1)...
          - exp(-omega*(N+Np)*delta)*fs_sequence_x(i+(N+Np)); %(1+exp(-omega))
      elseif Prova == 1
          beq_x = x(i) + xd(i)/omega - (1-exp(-omega*N*delta))*zx(i)+(1/omega)*integral...
              - omega*delta*exp(-omega*N*delta)*exp(-omega*delta*(0:Np-1))*fs_sequence_x(i+N:i+N+Np-1)...
              - exp(-omega*(N+Np)*delta)*fs_sequence_x(i+(N+Np)); %(1+exp(-omega))
      elseif Prova == 2
          beq_x = x(i) + xd(i)/omega - (1-exp(-omega*N*delta))*zx(i)+(1/omega)*integral+(1/omega)*knownDisturbance*(1/(1-exp(-omega*delta)))...
              - omega*delta*exp(-omega*N*delta)*exp(-omega*delta*(0:Np-1))*fs_sequence_x(i+N:i+N+Np-1)...
              - exp(-omega*(N+Np)*delta)*fs_sequence_x(i+(N+Np)); %(1+exp(-omega))
      elseif Prova == 3
          beq_x = x(i) + xd(i)/omega - (1-exp(-omega*N*delta))*zx(i)+(1/omega)*knownDisturbance*((1-exp(-omega*(N+1)*delta)))/(1-exp(-omega*delta))...
              - omega*delta*exp(-omega*N*delta)*exp(-omega*delta*(0:Np-1))*fs_sequence_x(i+N:i+N+Np-1)...
              - exp(-omega*(N+Np)*delta)*fs_sequence_x(i+(N+Np)); %(1+exp(-omega))
      elseif Prova == 4
          beq_x = x(i) + xd(i)/omega - (1-exp(-omega*N*delta))*zx(i)+(1/omega)*integral...
              - omega*delta*exp(-omega*N*delta)*exp(-omega*delta*(0:Np-1))*fs_sequence_x(i+N:i+N+Np-1)...
              - exp(-omega*(N+Np)*delta)*fs_sequence_x(i+(N+Np)); %(1+exp(-omega))
      elseif Prova == 5
          beq_x = x(i) + xd(i)/omega - (1-exp(-omega*N*delta))*zx(i)+(1/omega)*integral...
              - omega*delta*exp(-omega*N*delta)*exp(-omega*delta*(0:Np-1))*fs_sequence_x(i+N:i+N+Np-1)...
              - exp(-omega*(N+Np)*delta)*fs_sequence_x(i+(N+Np)); %(1+exp(-omega))
      end  
      
      %% QP solver
      zd_x = quadprog(H,f_x,A,b_x,Aeq,beq_x);
      
      z_pred_x = P*zd_x + zx(i);
      
      %% Model update
      x_updated = A_upd*[x(i); xd(i); zx(i)] + B_upd*(zd_x(1))+[0;delta;0]*F(i);
      
      x_obs(:,i+1) = x_obs(:,i) + delta*(A_obs*x_obs(:,i)+B_obs*zd_x(1)+L'*([x(i);zx(i)]-[x_obs(1,i);x_obs(3,i)]));  

      
      
      
      x(i+1) = x_updated(1);
      xd(i+1) = x_updated(2);
      zx(i+1) = x_updated(3);
      
      
      %% Divergent component of motion
      
      x_u(i) = x(1,i+1)+xd(1,i+1)/omega;
      
      %% Disturbance observer
      
      des_x(i+1)=((xd(i+1))-[0 1 0]*(A_upd*[x(i); xd(i); zx(i)] + B_upd*(zd_x(1))));
      dx = des_x(i);
      
      %% Feasibility region computation
      if Prova == 0
          x_u_min(i) = omega*Deltas*(zx_min(i:i+N-1)+restriction_x(1:N)) - (1/omega)*(delta*DDXX)*(1/(1-exp(-omega*delta))); % equivalent to -(1/omega)*A_disturbancex/delta (continous time!)
          x_u_max(i) = omega*Deltas*(zx_max(i:i+N-1)-restriction_x(1:N)) - (1/omega)*(delta*DDXX)*(1/(1-exp(-omega*delta)));
          x_u_min_nominal(i) = omega*Deltas*(zx_min(i:i+N-1));
          x_u_max_nominal(i) = omega*Deltas*(zx_max(i:i+N-1));
      elseif Prova == 1
          x_u_min(i) = omega*Deltas*(zx_min(i:i+N-1)) - (1/omega)*integral; % equivalent to -(1/omega)*A_disturbancex/delta (continous time!)
          x_u_max(i) = omega*Deltas*(zx_max(i:i+N-1)) - (1/omega)*integral;
          x_u_min_nominal(i) = omega*Deltas*(zx_min(i:i+N-1));
          x_u_max_nominal(i) = omega*Deltas*(zx_max(i:i+N-1));
      elseif Prova == 2
          x_u_min(i) = omega*Deltas*(zx_min(i:i+N-1)+restriction_x(1:N)) - 0*(1/omega)*integral- (1/omega)*knownDisturbance*(1/(1-exp(-omega*delta))); % equivalent to -(1/omega)*A_disturbancex/delta (continous time!)
          x_u_max(i) = omega*Deltas*(zx_max(i:i+N-1)-restriction_x(1:N)) - 0*(1/omega)*integral- (1/omega)*knownDisturbance*(1/(1-exp(-omega*delta)));
          x_u_min_nominal(i) = omega*Deltas*(zx_min(i:i+N-1));
          x_u_max_nominal(i) = omega*Deltas*(zx_max(i:i+N-1));
      elseif Prova == 3
          x_u_min(i) = omega*Deltas*(zx_min(i:i+N-1)) -  (1/omega)*knownDisturbance*((1-exp(-omega*(N+1)*delta)))/(1-exp(-omega*delta)); % equivalent to -(1/omega)*A_disturbancex/delta (continous time!)
          x_u_max(i) = omega*Deltas*(zx_max(i:i+N-1)) -  (1/omega)*knownDisturbance*((1-exp(-omega*(N+1)*delta)))/(1-exp(-omega*delta));
          x_u_min_nominal(i) = omega*Deltas*(zx_min(i:i+N-1));
          x_u_max_nominal(i) = omega*Deltas*(zx_max(i:i+N-1));
      elseif Prova == 5
          x_u_min(i) = omega*Deltas*(zx_min(i:i+N-1)+restriction_x(1:N)) - (1/omega)*integral; % equivalent to -(1/omega)*A_disturbancex/delta (continous time!)
          x_u_max(i) = omega*Deltas*(zx_max(i:i+N-1)-restriction_x(1:N)) - (1/omega)*integral;
          x_u_min_nominal(i) = omega*Deltas*(zx_min(i:i+N-1));
          x_u_max_nominal(i) = omega*Deltas*(zx_max(i:i+N-1));
                elseif Prova == 4
          x_u_min(i) = omega*Deltas*(zx_min(i:i+N-1)+restriction_x(1:N)) - (1/omega)*integral; % equivalent to -(1/omega)*A_disturbancex/delta (continous time!)
          x_u_max(i) = omega*Deltas*(zx_max(i:i+N-1)-restriction_x(1:N)) - (1/omega)*integral;
          x_u_min_nominal(i) = omega*Deltas*(zx_min(i:i+N-1));
          x_u_max_nominal(i) = omega*Deltas*(zx_max(i:i+N-1));
      end
      %        [zx_max, zx_min, fs_sequence_x] = FeasibilityFeedback(zx_max, zx_min, x_u(i), x_u_max(i), x_u_min(i), i, fs_sequence_x);

      
      i
     
      
      
%         figure(2)
%       clf
%       hold on
%       grid on
%       axis([0 1200 -0.5 0.5])
%       plot(disturbance(1:i+1)/delta, 'm', 'lineWidth', 2)
%       plot(x_obs(4,:), 'b', 'lineWidth', 2)
%       legend('Disturbance', 'Observed disturbance')
%       xlabel('time steps, \delta= 0.01 s')
%       ylabel('d [m/s^2]')
%       
%            drawnow
 
      
      
      
      
      
      end
      %% Plots
      
      clf 
            set(gca,'fontsize',12) 

      figure (1)
      clf
  
      hold on
      grid on
%       title('\fontsize{16}\color[rgb]{0 .5 .5}CoM and ZMP trajectories: C = 20, delta = 0.01') %: improving feasibility with step corrections
      plot(x,'r','lineWidth',2);
      plot(zx,'b','lineWidth',2);
%       plot(i:i+N-1,z_pred_x,'g','lineWidth',2);
%       plot(i+N:i+N+Np-1,fs_sequence_x(i+N:i+N+Np-1),'k','lineWidth',2);
      plot(zx_max,'m','lineWidth',2);
%       plot(fs_sequence_x,'y','lineWidth',1);
      plot(zx_min,'m','lineWidth',2);
      legend({'x','zx'},'FontSize',12)
      if balance == true
          axis([0 400 -2*(0.3*w+w) 2*(0.3*w+w)])
      else
          axis([0 400 -0.2 1])
      end
                             xlabel('time [s]')
                        ylabel('x [m]')
                        
      figure(3)
      clf
      grid on
      hold on
%       title('\fontsize{16}\color[rgb]{0 .5 .5}Divergent component of motion and feasibility region')
      plot(x_u,'m','lineWidth',2);
      plot(x_u_min,'b','lineWidth',2);
%       plot(x_u_min_nominal,'g')
%       plot(0.98*w*ones(1,i),'-.y')
%       plot(-0.98*w*ones(1,i),'-.y')
      plot(x_u_max,'b','lineWidth',2);
%       plot(x_u_max_nominal,'g')
%       axis([0 2000 -0.1 0.1])
      axis([0 400 -0.2 1])
%       set(gca,'XTickLabel',['0.0';'0.5';'1.0';'1.5';'2.0'; ...
%                             '2.5';'3.0';'3.5';'4.0'; '1.8'; ...
%                             '2.0';'2.2';'2.4';'2.6';'2.8'; ...
%                             '3.0';'3.2';'3.4';'3.6'; '3.8'; '4.0']) 
                        xlabel('time [s]')
                        ylabel('x_u [m]')

      legend({'x_u','feasibility region'},'FontSize',12)
      
      drawnow

      
%       figure(2)
%       clf
%       hold on
%       grid on
%       
%       for i=1:800
%       Dw_max(i) = w*(1-exp(-omega*i*delta))/(2*i*delta) - 0*exp((Np-i)*delta)*0.1/omega;
%       end
%       
%       plot(Dw_max,'r','lineWidth',2)
      
      figure(2)
      clf
      hold on
      grid on
      axis([0 400 -0.01 0.5])
      plot(disturbance, 'm', 'lineWidth', 2)
      plot(x_obs(4,:), 'b', 'lineWidth', 2)
      legend('Disturbance', 'Observed disturbance')
      xlabel('time [s]')
      ylabel('d [m/s^2]')
                 
      