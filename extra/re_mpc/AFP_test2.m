clear; clc; close all;

% Initialization and definition
% ---------------------------------------

delta = 0.01;
pre_time = 1;
h_com = 0.33;
g = 9.8;
omega = sqrt(g/h_com);
ss_time = 0.3;
ds_time = 0.2;
% w = 0.04;
w = 0.05;

first_support_foot = 'left'; %doesn't work

sim_time = 2.5;
N_sim = round(sim_time/delta);

%x, xd, xdd, y, yd, ydd, theta, time
initial_state = [0, 0, 0, 0, 0, 0, 0, 0];
observer_init = [0, 0, 0, 0, 0, 0, 0, 0];

solv = AFP_solver(initial_state, observer_init, pre_time, delta, omega, w, ss_time, ds_time, first_support_foot);
%solv.set_vref(0.1,0,(pi/16)/(ss_time+ds_time));
solv.set_vref(0.0,-1,0.0);

plot_options = PlotOptions();
plot_options.plot_com = 1;
plot_options.plot_zmp = 1;
plot_options.plot_pred_zmp = 0;
plot_options.plot_footsteps = 1;
plot_options.plot_pred_footsteps = 1;
plot_options.plot_pred_zmp_constraints = 1;
plot_options.plot_pred_footstep_constraints = 1;
plot_options.plot_orientation = 1;
solv.set_plot_options(plot_options);

% Set Gains: QZd, QJ, QVx, QVy, QZx, QZy, Qfsx, Qfsy, Q_omega, Q_obstacle
solv.set_gains(1,0,10,10,0,0,0,0,1,0);

solv.set_plot_limits(-0.1,0.5,-0.1,0.1);

% preassigned_footsteps_matrix = [    0, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04 ; ...
%                                  -0.1,  0.1, -0.1,  0.1, -0.1,  0.1, -0.1 ]';
% solv.set_footsteps( preassigned_footsteps_matrix );

% Simulation cycle
% ---------------------------------------

disp('Simulation cycle')

for i = 1 : N_sim 
    solv.cycle(i);
end




h = figure(1)
    
    % plot
    clf
    hold on
%     hold off
        axis([-0.05 0.4 -0.2 0.2])
        rect_x = [0.025,0.025,-0.025,-0.025,0.025];
        rect_y = [0.025,-0.025,-0.025,0.025,0.025];
        for k = 1:size(solv.footsteps,2)
        plot(solv.footsteps(1,k)+rect_x,solv.footsteps(2,k)+rect_y,'m','linewidth',2,'Handlevisibility','off')
        end
     solv.plot(i);  
     plot(solv.x_store,solv.y_store,'r','linewidth',2)
          plot(solv.zx_store,solv.zy_store,'b','linewidth',2)

        grid on
%                 axis([-0.05 0.3 -0.1 0.2])
             %0.05 0.15 -0.1 0.2])

      xlabel('x [m]')
      ylabel('y [m]')
      lgd =  legend('CoM', 'ZMP')
      lgd.FontSize = 12;
      drawnow

      %{
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'AFP_ReactiveStepping.pdf','-dpdf','-r0') 
%}



      
      set(gca,'fontsize',12)
      
      figure (2)
      clf
      hold on
      grid on
      plot(solv.x_store,'r','lineWidth',2);
      plot(solv.zx_store,'b','lineWidth',2);
      zx_min = [];
      zx_max = [];
      for k = 1:size(solv.footsteps,2)
          if k<size(solv.footsteps,2)
              zx_min = [zx_min,solv.footsteps(1,k)-ones(1,20)*w/2, (solv.footsteps(1,k)-w/2)+ (1:10)*(solv.footsteps(1,k+1)-solv.footsteps(1,k))/10];
              zx_max = [zx_max,solv.footsteps(1,k)+ones(1,20)*w/2, (solv.footsteps(1,k)+w/2)+ (1:10)*(solv.footsteps(1,k+1)-solv.footsteps(1,k))/10];
          else
              zx_min = [zx_min,solv.footsteps(1,k)-ones(1,20)*w/2];
              zx_max = [zx_max,solv.footsteps(1,k)+ones(1,20)*w/2];
          end
          
          %                 plot(solv.footsteps(1,k)-ones(20)*w/2,'m','linewidth',2,'Handlevisibility','off')
          %                 plot(solv.footsteps(1,k)+ones(20)*w/2,'m','linewidth',2,'Handlevisibility','off')
      end
      plot(zx_max,'m','lineWidth',2);
      plot(zx_min,'m','lineWidth',2);
      legend({'x','zx'},'FontSize',12)
      
      xlabel('time [s]')
      ylabel('x [m]')
      

      