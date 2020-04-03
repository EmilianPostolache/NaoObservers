classdef WallObstacle
    %OBSTACLE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        center
        clear_distance
        obstacle_is_active
    end
    
    methods
        function obj = WallObstacle(center, clear_distance)
            obj.center = center;
            obj.clear_distance = clear_distance;
            obj.obstacle_is_active = 1;
        end
        
        function versor = get_versor(obj, robot_pos)
            versor = [-1;0];
        end
        
        function cp = get_closest_point(obj, robot_pos)
            cp = [obj.center(1);robot_pos(2)];
        end
        
        function tp = get_tangent_point(obj, robot_pos)
            tp = obj.get_closest_point(robot_pos) + obj.get_versor(robot_pos)*obj.clear_distance;
        end
        
        function theta_normal = get_theta_normal(obj, robot_pos)
            theta_normal = -pi;
        end
        
        function theta_obstacle = get_theta_obstacle(obj, robot_pos)
            %DEPRECATED
            theta_obstacle = pi/2;
        end
        
        function distance = get_distance(obj, robot_pos)
            distance = abs(robot_pos(1)-obj.center(1));
        end
        
        function state = is_active(obj)
            if obj.obstacle_is_active == 1
                state = 1;
            else
                state = 0;
            end
        end
        
        function plot(obj, robot_pos)
            p=patch([obj.center(1) obj.center(1) 10 10],[-10 10 10 -10],'k');
            set(p,'FaceAlpha',0.1,'Linestyle','none');
            
            closest_point = obj.get_closest_point(robot_pos);
            tangent_point = obj.get_tangent_point(robot_pos);
            
            b=patch([tangent_point(1) tangent_point(1) obj.center(1) obj.center(1)],[-10 10 10 -10],[255 255 204]/255);
%             b=patch([tangent_point(1) tangent_point(1) 10 10],[-10 10 10 -10],[255 255 204]/255);
            set(b,'FaceAlpha',0.84,'Linestyle','none');
            
            if obj.is_active()
                %viscircles(closest_point',obj.clear_distance,'Color','r');
                line([tangent_point(1)-cos(obj.get_theta_obstacle(robot_pos)) tangent_point(1)+cos(obj.get_theta_obstacle(robot_pos))], ...
                     [tangent_point(2)-sin(obj.get_theta_obstacle(robot_pos)) tangent_point(2)+sin(obj.get_theta_obstacle(robot_pos))], ...
                     'color','k','LineWidth',1.5,'LineStyle',':')
                line([closest_point(1)-cos(obj.get_theta_obstacle(robot_pos)) closest_point(1)+cos(obj.get_theta_obstacle(robot_pos))], ...
                     [closest_point(2)-sin(obj.get_theta_obstacle(robot_pos)) closest_point(2)+sin(obj.get_theta_obstacle(robot_pos))], ...
                     'color','k','LineWidth',1.5,'LineStyle',':')
            end

%             plot(round_closest_point(1),round_closest_point(2),'*')
%             plot(round_tangent_point(1),round_tangent_point(2),'*')
        end
    end
    
end

