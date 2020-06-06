classdef PlotOptions
    %PLOT_OPTIONS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        plot_com
        plot_zmp
        plot_orientation
        plot_footsteps
        plot_pred_com
        plot_pred_zmp
        plot_pred_footsteps
        plot_pred_zmp_constraints
        plot_pred_footstep_constraints
    end
    
    methods
        function obj = PlotOptions()
            obj.plot_com = 1;
            obj.plot_zmp = 1;
            obj.plot_orientation = 1;
            obj.plot_footsteps = 1;
            obj.plot_pred_com = 1; %NOT AVAILABLE
            obj.plot_pred_zmp = 1;
            obj.plot_pred_footsteps = 1;
            obj.plot_pred_zmp_constraints = 1;
            obj.plot_pred_footstep_constraints = 1;
        end
    end
    
end

