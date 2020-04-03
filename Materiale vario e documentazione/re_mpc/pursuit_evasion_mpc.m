clear; clc; close all;

% Initialization and definition
% ---------------------------------------

delta = 0.05;
pre_time = 0.6;
h_com = 0.26;
g = 9.8;
omega = sqrt(g/h_com);
ss_time = 0.2;
ds_time = 0.1;
w = 0.04;
first_support_foot = 'right';

sim_time = 30;
N_sim = round(sim_time/delta);

%x, xd, xdd, y, yd, ydd, theta, time
initial_state_pursuer = [-0.4, 0, 0, 0.0, 0, 0, 0, 0];
initial_state_evader  = [ 0.13, 0, 0, 0.1, 0, 0, -pi, 0];

pursuer = AFP_solver(initial_state_pursuer, pre_time, delta, omega, w, ss_time, ds_time, first_support_foot);
evader  = AFP_solver(initial_state_evader , pre_time, delta, omega, w, ss_time, ds_time, first_support_foot);

% Set Gains: QZd, QJ, QVx, QVy, QZx, QZy, Qfsx, Qfsy, Q_omega, Q_obstacle
pursuer.set_gains(1,0,10,0,0,0,0,0,1,0.005);
evader.set_gains(1,0,10,0,0,0,0,0,1,0.005);

%gain 0.05 works fine with convex obstacles

% pursuer.set_plot_limits(-1.5, 1.5, -1.5, 1.5);
% evader.set_plot_limits(-1.5, 1.5, -1.5, 1.5);

% pursuer.set_plot_limits(-0.05, 0.45, -0.25, 0.25);
% evader.set_plot_limits(-0.05, 0.45, -0.25, 0.25);

plot_options = PlotOptions();
plot_options.plot_com = 1;
plot_options.plot_zmp = 1;
plot_options.plot_pred_zmp = 0;
plot_options.plot_footsteps = 1;
plot_options.plot_pred_footsteps = 1;
plot_options.plot_pred_zmp_constraints = 1;
plot_options.plot_pred_footstep_constraints = 1;
plot_options.plot_orientation = 1;
pursuer.set_plot_options(plot_options);
evader.set_plot_options(plot_options);

% Simulation cycle
% ---------------------------------------

disp('Simulation cycle')

% Gain
k = 0.8;

figure

for i = 1 : N_sim  
    % Compute theta_eva
    theta_aim = atan2(evader.y_store(i) - pursuer.y_store(i) , evader.x_store(i) - pursuer.x_store(i));
    pursuer_angle_rel_to_evader = wrapToPi(angdiff(theta_aim+pi,evader.theta));
    sign_evasion = -sign(pursuer_angle_rel_to_evader);
    theta_eva = atan2(pursuer.y_store(i) - evader.y_store(i) , pursuer.x_store(i) - evader.x_store(i)) + sign_evasion*pi/2;    
    
    % Set reference velocity
    pursuer.set_vref(0.1,0,k*angdiff(theta_aim,pursuer.theta));
%     pursuer.set_vref(0,0,0);
    evader.set_vref(-0.1,0,k*angdiff(theta_eva,evader.theta));
%     evader.set_vref(-0.1,0,0);
    
    % Obstacle parameters
    clear_distance = 0.07;
    pursuer_pos = [pursuer.x_store(i);pursuer.y_store(i)];
    evader_pos = [evader.x_store(i);evader.y_store(i)];
    
%     obstacle = [Obstacle([0.4; 0.7],0.03,clear_distance),...
%                 Obstacle([0.3; -0.3],0.02,clear_distance),
%                 WallObstacle([0.3; -0.0],clear_distance)];
            
    obstacle = [WallObstacle([0.4; -0.0],clear_distance)];
    
    closest_obstacle_distance_pursuer = 1000000;
    closest_obstacle_index_pursuer = 0;
    closest_obstacle_distance_evader = 1000000;
    closest_obstacle_index_evader = 0;
    
    for j = 1:size(obstacle,2)
        if obstacle(j).get_distance(pursuer_pos) < closest_obstacle_distance_pursuer
            closest_obstacle_index_pursuer = j;
            closest_obstacle_distance_pursuer = obstacle(j).get_distance(pursuer_pos);
        end
        if obstacle(j).get_distance(evader_pos) < closest_obstacle_distance_evader
            closest_obstacle_index_evader = j;
            closest_obstacle_distance_evader = obstacle(j).get_distance(evader_pos);
        end
        obstacle(j).obstacle_is_active = 0;
    end
    
    closest_point_evader = obstacle(closest_obstacle_index_evader).get_closest_point(evader_pos);
    tangent_point_evader = obstacle(closest_obstacle_index_evader).get_tangent_point(evader_pos);
    distance_evader = closest_obstacle_distance_evader;
    
    sign_test_evader = sign(wrapToPi(angdiff(obstacle(closest_obstacle_index_evader).get_theta_normal(evader_pos),evader.theta)));
    if sign_test_evader == 0
        sign_test_evader = 1;
    end
    theta_align_evader = wrapToPi(obstacle(closest_obstacle_index_evader).get_theta_normal(evader_pos) - sign_test_evader*pi/2);
    
    evader.add_obstacle_constraint_on_footsteps(closest_point_evader, tangent_point_evader ,theta_align_evader);
    obstacle(closest_obstacle_index_evader).obstacle_is_active = 1;
    
    closest_point_pursuer = obstacle(closest_obstacle_index_pursuer).get_closest_point(pursuer_pos);
    tangent_point_pursuer = obstacle(closest_obstacle_index_pursuer).get_tangent_point(pursuer_pos);
    distance_pursuer = closest_obstacle_distance_pursuer;
    
    sign_test_pursuer = sign(wrapToPi(angdiff(obstacle(closest_obstacle_index_pursuer).get_theta_normal(pursuer_pos),pursuer.theta)));
    if sign_test_pursuer == 0
        sign_test_pursuer = 1;
    end
    theta_align_pursuer = wrapToPi(obstacle(closest_obstacle_index_pursuer).get_theta_normal(pursuer_pos) - sign_test_pursuer*pi/2);
    
    pursuer.add_obstacle_constraint_on_footsteps(closest_point_pursuer, tangent_point_pursuer ,theta_align_pursuer);
    obstacle(closest_obstacle_index_pursuer).obstacle_is_active = 1;
    
    % Cycle
    pursuer.cycle(i);
    evader.cycle(i);
    
    % Plots
    clf
    hold on
    quiver(pursuer_pos(1),pursuer_pos(2),0.1*cos(theta_align_pursuer),0.1*sin(theta_align_pursuer),'LineWidth',1,'MaxHeadSize',2,'color','g')
    quiver(evader_pos(1),evader_pos(2),0.1*cos(theta_align_evader),0.1*sin(theta_align_evader),'LineWidth',1,'MaxHeadSize',2,'color','g')

    
    % Plot obstacles
    for j = 1:size(obstacle,2)
%         obstacle(j).plot(pursuer_pos);
        obstacle(j).plot(evader_pos);
    end
       
    % Plot the robots
    pursuer.plot(i);
    evader.plot(i);

    hold off
    grid on
    
    axis([-0.5,0.5,-0.3,0.3])
    
    drawnow
end