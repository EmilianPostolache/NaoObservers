classdef Obstacle
    %OBSTACLE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        center
        size
        clear_distance
        obstacle_is_active
    end
    
    methods
        function obj = Obstacle(center, size, clear_distance)
            obj.center = center;
            obj.size = size;
            obj.clear_distance = clear_distance;
            obj.obstacle_is_active = 1;
        end
        
        function versor = get_versor(obj, robot_pos)
            versor = (robot_pos-obj.center) / sqrt((robot_pos(1)-obj.center(1))^2+(robot_pos(2)-obj.center(2))^2);
        end
        
        function cp = get_closest_point(obj, robot_pos)
            versor = (robot_pos-obj.center) / sqrt((robot_pos(1)-obj.center(1))^2+(robot_pos(2)-obj.center(2))^2);
            cp = obj.center + versor*obj.size;
        end
        
        function tp = get_tangent_point(obj, robot_pos)
            versor = (robot_pos-obj.center) / sqrt((robot_pos(1)-obj.center(1))^2+(robot_pos(2)-obj.center(2))^2);
            tp = obj.center + versor*(obj.size+obj.clear_distance);
        end
        
        function theta_normal = get_theta_normal(obj, robot_pos)
            normal_versor = -(robot_pos-obj.center) / sqrt((robot_pos(1)-obj.center(1))^2+(robot_pos(2)-obj.center(2))^2);
            theta_normal = atan2(normal_versor(2),normal_versor(1));
        end
        
        function theta_obstacle = get_theta_obstacle(obj, robot_pos)
            versor = (robot_pos-obj.center) / sqrt((robot_pos(1)-obj.center(1))^2+(robot_pos(2)-obj.center(2))^2);
            tangent_versor = [0,1;-1,0]*versor;
            theta_obstacle = atan2(tangent_versor(2),tangent_versor(1));
        end
        
        function distance = get_distance(obj, robot_pos)
            tangent_point = obj.get_tangent_point(robot_pos);
            distance = sqrt((robot_pos(1)-tangent_point(1))^2 + (robot_pos(2)-tangent_point(2))^2);
        end
        
        function state = is_active(obj)
            if obj.obstacle_is_active == 1
                state = 1;
            else
                state = 0;
            end
        end
        
        function plot(obj, robot_pos)
            closest_point = obj.get_closest_point(robot_pos);
            tangent_point = obj.get_tangent_point(robot_pos);
            
            patch([obj.center(1)+obj.size obj.center(1)+obj.size ...
                   obj.center(1)-obj.size obj.center(1)-obj.size], ...
                  [obj.center(2)+obj.size obj.center(2)-obj.size ...
                   obj.center(2)-obj.size obj.center(2)+obj.size],'b');
            if obj.is_active()
%                 viscircles(closest_point',obj.clear_distance,'Color','r');
                line([tangent_point(1)-cos(obj.get_theta_obstacle(robot_pos)) tangent_point(1)+cos(obj.get_theta_obstacle(robot_pos))], ...
                     [tangent_point(2)-sin(obj.get_theta_obstacle(robot_pos)) tangent_point(2)+sin(obj.get_theta_obstacle(robot_pos))])
            end

%             plot(round_closest_point(1),round_closest_point(2),'*')
%             plot(round_tangent_point(1),round_tangent_point(2),'*')
        end
    end
    
end

