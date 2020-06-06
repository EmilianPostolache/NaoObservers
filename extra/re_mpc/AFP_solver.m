classdef AFP_solver < dynamicprops
    
    properties(GetAccess=public)
        % gains
        QZd
        QJ
        QVx
        QVy
        QZx
        QZy
        Qfsx
        Qfsy
        Q_omega
        Q_obstacle
        
        % parameters
        ss_time
        ds_time
        S_samples
        D_samples
        M_samples
        N_samples
        F_samples
        delta
        omega
        w
        preassigned_footsteps_matrix
        footsteps_are_preassigned
        first_support_foot
        box_a
        box_b
        box_c
        
        % state
        x
        xd
        y
        yd
        t_curr
        zx
        zy
        theta
        iter
        iterF
        fxc
        fyc
        footstep
        vref_x
        vref_y
        vref_omega
        x_hat
        xd_hat
        zx_hat
        d_x_hat
        y_hat
        yd_hat
        zy_hat
        d_y_hat
        xu_Max
        xu_min
        yu_Max
        yu_min
        
        % plot utilities
        zd
        pred_zx
        pred_zy
        pred_fs
        plot_limits
        plot_options
        
        % storage arrays
        x_store
        y_store
        xd_store
        yd_store
        zx_store
        zy_store
        zxd_store
        zyd_store
        footsteps
        
        % solving utilities
        HQ
        fQ
        Vu
        Vs
        Cc
        A_zmp
        b_zmp
        A_fs
        b_fs
        Aeq
        beq
        predicted_rotations
        
        % obstacles
        A_obstacle
        b_obstacle
        obstacle_is_on
        obstacle_theta_abs
        obstacle_theta
        obstacle_closest_point_abs
        obstacle_closest_point
        obstacle_tangent_point_abs
        obstacle_tangent_point
        
        % deltas
        DeltasN
        DeltasP
        
        
        
    end
    
    methods
        function obj = AFP_solver(state, observer_init, pre_time, delta, omega, w, ss_time, ds_time, first_support_foot)
            
            % Initialize global variables
            obj.delta = delta;
            obj.omega = omega;
            obj.w = w;
            obj.ss_time = ss_time;
            obj.ds_time = ds_time;
            
            % Initialize state
            state = num2cell(state');
            [obj.x, obj.xd, xdd, obj.y, obj.yd, ydd, obj.theta, obj.t_curr] = state{:};
            state_obs = num2cell(observer_init');
            [obj.x_hat, obj.xd_hat, obj.zx_hat, obj.d_x_hat, obj.y_hat, obj.yd_hat, obj.zy_hat, obj.d_y_hat] = state_obs{:};
            
            
            first_foot_position = [cos(obj.theta),-sin(obj.theta);sin(obj.theta),cos(obj.theta)]*[0;0.05];
            %first_foot_position = [0;0.05];
            
            obj.fxc = obj.x + first_foot_position(1);
            obj.fyc = obj.y + first_foot_position(2);
            
            %obj.x = obj.x - obj.fxc;
            %obj.y = obj.y - obj.fyc;
            
            obj.x = 0;
            obj.y = -0.05;
            
            obj.zx = obj.x - (1/obj.omega^2)*xdd;
            obj.zy = obj.y - (1/obj.omega^2)*ydd;
            
            obj.N_samples = round(pre_time/obj.delta);
            obj.S_samples = round(obj.ss_time/obj.delta);
            obj.D_samples = round(obj.ds_time/obj.delta);
            obj.F_samples = (obj.S_samples+obj.D_samples)*2;
            obj.M_samples = ceil(obj.N_samples/(obj.S_samples+obj.D_samples));
            N=obj.N_samples;
            S=obj.S_samples;
            D=obj.D_samples;
            F=obj.F_samples;
            M=obj.M_samples;
            
            % useful vector to compute feasibility region
            
            for j = 1:N
                
                obj.DeltasN(j) = obj.omega*obj.delta*exp(-omega*delta*(j-1));
                obj.DeltasP(j) = obj.omega*obj.delta*exp(-omega*delta*(j+N-1));
                    
            end

                        
            % Initialize storage vectors           
            rot = [cos(obj.theta),-sin(obj.theta);sin(obj.theta),cos(obj.theta)];
            
            pos_abs = [obj.fxc;obj.fyc] + rot*[obj.x;obj.y];
            vel_abs = rot*[obj.xd;obj.yd];
            zmp_abs = [obj.fxc;obj.fyc] + rot*[obj.zx;obj.zy];
            obj.x_store(1) = pos_abs(1);
            obj.y_store(1) = pos_abs(2);
            obj.xd_store(1) = vel_abs(1);
            obj.yd_store(1) = vel_abs(2);
            obj.zx_store(1) = zmp_abs(1);
            obj.zy_store(1) = zmp_abs(2);
            
            % Velocity tracking matrices
            
            ch = cosh(obj.omega*obj.delta);
            sh = sinh(obj.omega*obj.delta);
            A_upd = [ch, sh/obj.omega, 1-ch; obj.omega*sh, ch, -obj.omega*sh; 0, 0, 1];
            B_upd = [obj.delta-sh/obj.omega; 1-ch; obj.delta];
            Vu = [];
            Vs = [];
            for i = 1:N
                Vs_newline = [obj.omega*sh, ch, -obj.omega*sh] * A_upd^(i-1);
                Vu_newline = zeros(1,N);
                Vu_newline(1,i) = 1-ch;
                if i > 1
                    for j = 0:(i-2)
                        Vu_newline(1,j+1) = [obj.omega*sh, ch, -obj.omega*sh] * A_upd^(i-j-2)*B_upd;
                    end
                end
                Vu = [Vu; Vu_newline];
                Vs = [Vs; Vs_newline];
            end
            
            obj.Vu = Vu;
            obj.Vs = Vs;
            
            obj.QZd = 1;
            obj.QVx = 0;
            obj.QVy = 0;
            obj.QZx = 0;
            obj.QZy = 0;
              
            obj.footstep = 0;
            obj.first_support_foot = first_support_foot;
            
            obj.vref_x = 0;
            obj.vref_y = 0;
            
            obj.box_a = 0.1;
            obj.box_b = 0.15; %0.15
            obj.box_c = 0.05;
            
            obj.plot_limits = [-1, 1, -1, 1];
            
            obj.footsteps(:,1) = [obj.fxc;obj.fyc];
            
            obj.footsteps_are_preassigned = 0;
            obj.preassigned_footsteps_matrix = [];
            
            obj.A_obstacle = [];
            obj.b_obstacle = [];
            
            obj.obstacle_is_on = 0;
            obj.obstacle_theta = 0;            
        end
        
        function set_footsteps(obj, preassigned_footsteps_matrix)
            obj.footsteps_are_preassigned = 1;
            obj.preassigned_footsteps_matrix = preassigned_footsteps_matrix;
        end
        
        function set_vref(obj, vx, vy, ang_vel)
            obj.vref_x = vx;
            obj.vref_y = vy;
            obj.vref_omega = ang_vel;
        end
        
        function set_gains(obj, QZd, QJ, QVx, QVy, QZx, QZy, Qfsx, Qfsy, Q_omega, Q_obstacle)
            obj.QZd = QZd;
            obj.QJ  = QJ;
            obj.QVx = QVx;
            obj.QVy = QVy;
            obj.QZx = QZx;
            obj.QZy = QZy;
            obj.Qfsx = Qfsx;
            obj.Qfsy = Qfsy;
            obj.Q_omega = Q_omega;
            obj.Q_obstacle = Q_obstacle;
        end

        function gen_zmp_constraints(obj)
            N=obj.N_samples;S=obj.S_samples;D=obj.D_samples;F=obj.F_samples;M=obj.M_samples;
            
            
            p = ones(N,1);
            P = tril(ones(N)) * obj.delta;
            
            Cstart = zeros(N,1);
            if obj.footstep == 0 && obj.iterF <= S
                Cstart(1:S-obj.iterF) = ones(S-obj.iterF,1);
            end
            
            %step_theta = obj.vref_omega*(obj.ss_time+obj.ds_time);
            
            Icf = blkdiag(eye(S),zeros(D,D));

            if obj.footstep == 0
              %  Ic = zeros(S+D-obj.iterF, S+D-obj.iterF);
            else
              %  Ic = Icf(obj.iterF+1:end,obj.iterF+1:end);
            end
            
            
            % ATTENTION : this should be the only part with S,D
            
            rot_cos = cos(0)*eye(F/2-obj.iterF-1);
            rot_sin = sin(0)*eye(F/2-obj.iterF-1);
            fullCc_step = [ones(S,1);fliplr(0:D-1)'/D];

            fullCc_nextstep = [zeros(S,1);(1:D)'/D];

%             Cc = zeros(F/2-obj.iterF,1);
            fullCc = [fullCc_nextstep, zeros(F/2,1); fullCc_step, fullCc_nextstep; zeros(F/2,1), fullCc_step];

            Cc = fullCc((obj.iterF+1):(obj.iterF + N), :);
            
            
            
            for i = 1:M-1
              %  Ic = blkdiag(Ic, Icf);
                rot_cos = blkdiag(rot_cos, cos(obj.predicted_rotations(i+1))*eye(F/2));
                rot_sin = blkdiag(rot_sin, sin(obj.predicted_rotations(i+1))*eye(F/2));
%                 Cc = blkdiag(Cc, ones(F/2,1));
%                 Cc = blkdiag(Cc, )
            end
            
           % Ic = blkdiag(Ic, Icf(1:obj.iterF,1:obj.iterF));
            rot_cos = blkdiag(rot_cos, cos(obj.predicted_rotations(M+1))*eye(obj.iterF+1));
            rot_sin = blkdiag(rot_sin, sin(obj.predicted_rotations(M+1))*eye(obj.iterF+1));
%             Cc = blkdiag(Cc, ones(obj.iterF,1));
                      
%             Cc = Cc(:,2:end);
            obj.Cc = Cc;
            
            % REMOVING Ic
            Ic = eye(N);


            rot = [ rot_cos,-rot_sin; ...
                    rot_sin, rot_cos  ];
            
            AQ_max =  rot'*[   Ic*P ,  -Ic*Cc , 0*Ic*P , -0*Ic*Cc; ...
                             0*Ic*P , 0*Ic*Cc ,   Ic*P ,   -Ic*Cc  ];
                        
            AQ_min = -rot'*[   Ic*P ,  -Ic*Cc , 0*Ic*P , -0*Ic*Cc; ...
                             0*Ic*P , 0*Ic*Cc ,   Ic*P ,   -Ic*Cc  ];
            
                         %FIX FIRST STEP
            first_step_constraint_max = [0*Ic*p;0.1*Ic*Cstart];
            first_step_constraint_min = [0*Ic*p;0.1*Ic*Cstart];
            
            bQ_max = [Ic*p;Ic*p]*obj.w/2 - rot'*[Ic*p*obj.zx;Ic*p*obj.zy] + first_step_constraint_max;
            bQ_min = [Ic*p;Ic*p]*obj.w/2 + rot'*[Ic*p*obj.zx;Ic*p*obj.zy] + first_step_constraint_min;
            
            obj.A_zmp = [AQ_max;AQ_min];
            obj.b_zmp = [bQ_max;bQ_min];
        end

        function gen_footstep_constraints(obj)
            N=obj.N_samples;S=obj.S_samples;D=obj.D_samples;F=obj.F_samples;M=obj.M_samples;
            
            pfR = zeros(M,1);
            pfL = zeros(M,1);
            if mod(obj.footstep,2) == 0
                for i = 1:M
                    if mod(i,2) == 0
                        pfR(i) = 1;
                    else
                        pfL(i) = 1;
                    end
                end
            else
                for i = 1:M
                    if mod(i,2) == 0
                        pfL(i) = 1;
                    else
                        pfR(i) = 1;
                    end
                end
            end
            
            rot_cos = zeros(M,M);
            rot_sin = zeros(M,M);
            %step_theta = obj.vref_omega*(obj.ss_time+obj.ds_time);

            rot_cos(1,1) = cos(0);
            rot_sin(1,1) = sin(0);
            
            for j = 2:M
                % Fill the rotation matrix
                rot_cos(j,j) = cos(obj.predicted_rotations(i));
                rot_sin(j,j) = sin(obj.predicted_rotations(i));
            end

            rot = [ rot_cos,-rot_sin; ...
                    rot_sin, rot_cos  ];
                
            A = eye(M) - [zeros(1,M);eye(M-1),zeros(M-1,1)];
            
            A_fs_max = rot'*[zeros(M,N) , A , zeros(M,N) , zeros(M,M); ...
                             zeros(M,N) , zeros(M,M) , zeros(M,N) , A]; ...
            A_fs_min = rot'*[zeros(M,N) ,-A , zeros(M,N) , zeros(M,M); ...
                             zeros(M,N) , zeros(M,M) , zeros(M,N) ,-A];
            
            obj.A_fs = [A_fs_max;A_fs_min];

            b_fs_max = [ones(M,1)*obj.box_c ; -pfL*obj.box_a + pfR*obj.box_b];
            b_fs_min = [ones(M,1)*obj.box_c ;  pfL*obj.box_b - pfR*obj.box_a];
            
            obj.b_fs = [b_fs_max;b_fs_min];
        end
        
        function gen_stability_constraint(obj)
            
            
            BD = 0;
            
            N=obj.N_samples;S=obj.S_samples;D=obj.D_samples;F=obj.F_samples;M=obj.M_samples;
            
            if N == F
                lambda = exp(-obj.omega*obj.delta);
                b_exp = exp(-obj.omega*obj.delta*(0:N-1))';

                obj.Aeq = [(1/obj.omega)*(1-lambda)/(1-lambda^N)*b_exp', zeros(1,M), zeros(1,N), zeros(1,M); ...
                           zeros(1,N), zeros(1,M) , (1/obj.omega)*(1-lambda)/(1-lambda^N)*b_exp', zeros(1,M)];

                obj.beq = [obj.x + obj.xd/obj.omega - obj.zx + BD*(1/obj.omega^2)*obj.d_x_hat; ...
                           obj.y + obj.yd/obj.omega - obj.zy + BD*(1/obj.omega^2)*obj.d_y_hat];
                   
            elseif N < F
                size_z = size(obj.zxd_store,2);
                
                linear_completion_on = 1;
                
                if (linear_completion_on && obj.iter > F-N) || not(linear_completion_on)
                    
                    if obj.iter > F-N
                        last_samples_x = obj.zxd_store(size_z-(F-N)+1:size_z)';
                        last_samples_y = obj.zyd_store(size_z-(F-N)+1:size_z)';
                    else
                        last_samples_x = zeros(F-N,1);
                        last_samples_y = zeros(F-N,1);
                    end
                    
                    L = exp(-obj.omega*obj.delta);
                    b_exp = exp(-obj.omega*obj.delta*(0:N-1))';
                    
                    obj.Aeq = [(1/obj.omega)*(1-L)/(1-L^F)*b_exp', zeros(1,M), zeros(1,N), zeros(1,M); ...
                        zeros(1,N), zeros(1,M) , (1/obj.omega)*(1-L)/(1-L^F)*b_exp', zeros(1,M)];
                    
                    obj.beq = [obj.x + obj.xd/obj.omega - obj.zx ...
                        - L^N*(1/obj.omega)*(1-L)/(1-L^F)*exp(-obj.omega*obj.delta*(0:F-N-1))*last_samples_x ; ...
                        obj.y + obj.yd/obj.omega - obj.zy ...
                        - L^N*(1/obj.omega)*(1-L)/(1-L^F)*exp(-obj.omega*obj.delta*(0:F-N-1))*last_samples_y];

                else
                    disp('Linear completion')
                    L = exp(-obj.omega*obj.delta);
                    b_exp = exp(-obj.omega*obj.delta*(0:N-1))';
                    coeff = ((F-N-1)*L^(F-N+1) - (F-N)*L*(F-N) + L)/((F-N)*(L-1)^2);
                    
                    add1 = [zeros(1,N-1),1]*(L^N*(1-L^(F-N)-(1-L)*coeff))/(obj.omega*(1-L^N));
                    add2 = [1,zeros(1,N-1)]*(L^N*(1-L)*coeff)/(obj.omega*(1-L^N));
                    
                    obj.Aeq = [(1/obj.omega)*(1-L)/(1-L^F)*b_exp'+add1+add2, zeros(1,M), zeros(1,N), zeros(1,M); ...
                        zeros(1,N), zeros(1,M) , (1/obj.omega)*(1-L)/(1-L^F)*b_exp'+add1+add2, zeros(1,M)];
                    
                    obj.beq = [obj.x + obj.xd/obj.omega - obj.zx; ...
                        obj.y + obj.yd/obj.omega - obj.zy ];
                end
            elseif N > F
                lambda = exp(-obj.omega*obj.delta);
                b_start  = [ exp(-obj.omega*obj.delta*(0:N-F-1))' ; zeros(F,1)       ];
                b_period = [ zeros(N-F,1) ; exp(-obj.omega*obj.delta*(0:F-1))'       ];

                A_start = [((1-lambda)/obj.omega)*b_start', zeros(1,M), zeros(1,N), zeros(1,M); ...
                           zeros(1,N), zeros(1,M) , ((1-lambda)/obj.omega)*b_start', zeros(1,M)];

                A_period = [(lambda^(N-F)/obj.omega)*(1-lambda)/(1-lambda^F)*b_period', zeros(1,M), zeros(1,N), zeros(1,M); ...
                            zeros(1,N), zeros(1,M) , (lambda^(N-F)/obj.omega)*(1-lambda)/(1-lambda^F)*b_period', zeros(1,M)];

                obj.Aeq = A_start + A_period;
                   
                obj.beq = [obj.x + obj.xd/obj.omega - obj.zx; ...
                           obj.y + obj.yd/obj.omega - obj.zy ];
            end
        end
        
        function gen_stability_constraint_constant(obj)
            N=obj.N_samples;S=obj.S_samples;D=obj.D_samples;F=obj.F_samples;M=obj.M_samples;
            
            lambda = exp(-obj.omega*obj.delta);
            b_exp = exp(-obj.omega*obj.delta*(0:N-1))';
            
            obj.Aeq = [((1-lambda)/obj.omega)*b_exp', zeros(1,M), zeros(1,N), zeros(1,M); ...
                      zeros(1,N), zeros(1,M) , ((1-lambda)/obj.omega)*b_exp', zeros(1,M)];
            
            obj.beq = [obj.x + obj.xd/obj.omega - obj.zx; ...
                       obj.y + obj.yd/obj.omega - obj.zy ];
        end
        
        function gen_stability_constraint_footsteps(obj)
            N=obj.N_samples;S=obj.S_samples;D=obj.D_samples;F=obj.F_samples;M=obj.M_samples;
            
            footstep_velocity_matrix = obj.preassigned_footsteps_matrix/obj.ds_time;
            preassigned_footsteps_number = size(footstep_velocity_matrix,1);
            
            lambda = exp(-obj.omega*obj.delta);
            b_exp = exp(-obj.omega*obj.delta*(0:N-1))';
            b_exp_fs = exp(-obj.omega*obj.delta*(S+D)*(0:preassigned_footsteps_number-1))';
            
            %prediction component
            obj.Aeq = [((1-lambda)/obj.omega)*b_exp' - lambda^N*ones(1,N)*obj.delta, zeros(1,M), zeros(1,N), zeros(1,M); ...
                      zeros(1,N), zeros(1,M) , ((1-lambda)/obj.omega)*b_exp' - lambda^N*ones(1,N)*obj.delta, zeros(1,M)];
            
            %footstep component
            footsteps_remaining = size(obj.preassigned_footsteps_matrix,1) - obj.footstep;
            predictable_footsteps = min(M,footsteps_remaining);
            
            %starting position of the tail
            fs_component_start_x = sum(footstep_velocity_matrix(1:predictable_footsteps,1))*obj.ds_time;
            fs_component_start_y = sum(footstep_velocity_matrix(1:predictable_footsteps,2))*obj.ds_time;
            
            if mod(obj.iterF+N, S+D) < S
                %prediction ends in single support
            else
                %prediction ends in double support
            end
            
            footstep_velocity_matrix(1:predictable_footsteps,:) = zeros(predictable_footsteps,2);

            obj.beq = [obj.x + obj.xd/obj.omega - obj.zx;  ...
                       obj.y + obj.yd/obj.omega - obj.zy ] ...
                      - lambda^(M*(S+D)+S-obj.iterF)*[((1-lambda^D)/obj.omega)*b_exp_fs'*footstep_velocity_matrix(:,1); ... 
                                                      ((1-lambda^D)/obj.omega)*b_exp_fs'*footstep_velocity_matrix(:,2)] ...
                      - lambda^N*[fs_component_start_x;fs_component_start_y];

        end

        function add_obstacle_constraint_on_footsteps(obj, closest_point, tangent_point, obstacle_theta)                                      
            obj.obstacle_closest_point_abs = closest_point;
            obj.obstacle_tangent_point_abs = tangent_point;
            obj.obstacle_theta_abs = obstacle_theta;
            obj.obstacle_is_on = 1;
        end
        
        function gen_obstacle_constraint_on_footsteps(obj)
            N=obj.N_samples;S=obj.S_samples;D=obj.D_samples;F=obj.F_samples;M=obj.M_samples;
            
            p = obj.obstacle_tangent_point;
            
            if 0 < -sin(obj.obstacle_theta)*(p(1)) + cos(obj.obstacle_theta)*(p(2))
                sign =  1;
            else
                sign = -1;
            end
                  
            obj.A_obstacle = [obj.A_obstacle; ...
                              sign*[zeros(M,N), -eye(M)*sin(obj.obstacle_theta), ...
                                    zeros(M,N),  eye(M)*cos(obj.obstacle_theta)]];

            obj.b_obstacle = [obj.b_obstacle; ...
                              sign*ones(M,1)*(-sin(obj.obstacle_theta)*(p(1))+ ...
                                               cos(obj.obstacle_theta)*(p(2)))];
                                           
        end

        function gen_cost_function(obj)
%             if obj.obstacle_is_on
%                 obj.HQ = obj.H_Zd() + obj.H_Xj() + obj.H_V() + obj.H_Z() + obj.H_footsteps() + obj.H_obstacle();
%                 obj.fQ = obj.f_Zd() + obj.f_Xj() + obj.f_V() + obj.f_Z() + obj.f_footsteps() + obj.f_obstacle();
%             else
                obj.HQ = obj.H_Zd() + obj.H_Xj() + obj.H_V() + obj.H_Z() + obj.H_footsteps();
                obj.fQ = obj.f_Zd() + obj.f_Xj() + obj.f_V() + obj.f_Z() + obj.f_footsteps();
%             end
        end
        
        function H = H_Zd(obj)
            N=obj.N_samples;S=obj.S_samples;D=obj.D_samples;F=obj.F_samples;M=obj.M_samples;

            HQ1 = zeros(N+M, N+M);
            HQ1(1:N, 1:N) = obj.QZd*eye(N);
            HQ2 = zeros(N+M, N+M);
            HQ2(1:N, 1:N) = obj.QZd*eye(N);
            H = blkdiag(HQ1,HQ2);
        end
        
        function f = f_Zd(obj)
            N=obj.N_samples;S=obj.S_samples;D=obj.D_samples;F=obj.F_samples;M=obj.M_samples;
            
            fQ1 = zeros(N+M,1);
            fQ2 = zeros(N+M,1);
            f = [fQ1;fQ2];
        end
        
        function H = H_Xj(obj)
            N=obj.N_samples;S=obj.S_samples;D=obj.D_samples;F=obj.F_samples;M=obj.M_samples;

            HQ1 = zeros(N+M, N+M);
            HQ1(1:N, 1:N) = obj.QJ*obj.omega^4 * ((obj.Vu-eye(N))'*(obj.Vu-eye(N)));
            HQ2 = zeros(N+M, N+M);
            HQ2(1:N, 1:N) = obj.QJ*obj.omega^4 * ((obj.Vu-eye(N))'*(obj.Vu-eye(N)));
            H = blkdiag(HQ1,HQ2);
        end
        
        function f = f_Xj(obj)
            N=obj.N_samples;S=obj.S_samples;D=obj.D_samples;F=obj.F_samples;M=obj.M_samples;
            
            fQ1 = obj.QJ*obj.omega^4 * [(obj.Vu-eye(N))'*(obj.Vs*[obj.x;obj.xd;obj.zx]); zeros(M,1)];
            fQ2 = obj.QJ*obj.omega^4 * [(obj.Vu-eye(N))'*(obj.Vs*[obj.y;obj.yd;obj.zy]); zeros(M,1)];
            f = [fQ1;fQ2];
        end
        
        function H = H_V(obj)
            N=obj.N_samples;S=obj.S_samples;D=obj.D_samples;F=obj.F_samples;M=obj.M_samples;

            HQ1 = zeros(N+M, N+M);
            HQ1(1:N, 1:N) = obj.QVx * (obj.Vu'*obj.Vu);
            HQ2 = zeros(N+M, N+M);
            HQ2(1:N, 1:N) = obj.QVy * (obj.Vu'*obj.Vu);
            H = blkdiag(HQ1,HQ2);
        end
        
        function f = f_V(obj)
            N=obj.N_samples;S=obj.S_samples;D=obj.D_samples;F=obj.F_samples;M=obj.M_samples;
            
            % Velocity arc to track
            for i = 1:N
                v_arc_x(i) =  obj.vref_x*cos(i*obj.vref_omega*obj.delta);
                v_arc_y(i) =  obj.vref_x*sin(i*obj.vref_omega*obj.delta);
            end
            
            fQ1 = obj.QVx*[obj.Vu'*(obj.Vs*[obj.x;obj.xd;obj.zx]-v_arc_x'); zeros(M,1)];
            fQ2 = obj.QVy*[obj.Vu'*(obj.Vs*[obj.y;obj.yd;obj.zy]-v_arc_y'); zeros(M,1)];
            f = [fQ1;fQ2];
        end
        
        function H = H_Z(obj)
            N=obj.N_samples;S=obj.S_samples;D=obj.D_samples;F=obj.F_samples;M=obj.M_samples;
            Cc = obj.Cc;
            p = ones(N,1);
            P = tril(ones(N)) * obj.delta;
            
            HQ1 = zeros(N+M, N+M);
            HQ1(1:N, 1:N) = obj.QZx * (P'*P);
            HQ1(1:N, N+1:N+M) = -obj.QZx * (P'*Cc);
            HQ1(N+1:N+M, 1:N) = -obj.QZx * (P'*Cc)';%(Cc'*P);
            HQ1(N+1:N+M, N+1:N+M) = obj.QZx * (Cc'*Cc);
            HQ2 = zeros(N+M, N+M);
            HQ2(1:N, 1:N) = obj.QZy * (P'*P);
            HQ2(1:N, N+1:N+M) = -obj.QZy * (P'*Cc);
            HQ2(N+1:N+M, 1:N) = -obj.QZy * (P'*Cc)';%(Cc'*P);
            HQ2(N+1:N+M, N+1:N+M) = obj.QZy * (Cc'*Cc);
            
            H = blkdiag(HQ1,HQ2);
        end

        function f = f_Z(obj)
            N=obj.N_samples;S=obj.S_samples;D=obj.D_samples;F=obj.F_samples;M=obj.M_samples;
            Cc = obj.Cc;
            p = ones(N,1);
            P = tril(ones(N)) * obj.delta;
            
            fQ1 = [ obj.QZx * P'  * p *obj.zx; ...
                   -obj.QZx * Cc' * p *obj.zx];
            fQ2 = [ obj.QZy * P'  * p *obj.zy; ...
                   -obj.QZy * Cc' * p *obj.zy];
              
            f = [fQ1;fQ2];
        end
        
        function H = H_footsteps(obj)
            N=obj.N_samples;S=obj.S_samples;D=obj.D_samples;F=obj.F_samples;M=obj.M_samples;

            I_diff = eye(M) - [zeros(1,M);eye(M-1),zeros(M-1,1)];
            
            H_fs1 = obj.Qfsx*(I_diff'*I_diff);
            H_fs2 = obj.Qfsy*(I_diff'*I_diff);
            
            H = [zeros(N,N),zeros(N,M),zeros(N,N),zeros(N,M); ...
                 zeros(M,N),  H_fs1   ,zeros(M,N),zeros(M,M); ...
                 zeros(N,N),zeros(N,M),zeros(N,N),zeros(N,M); ...
                 zeros(M,N),zeros(M,M),zeros(M,N),  H_fs2  ];
        end
        
        function f = f_footsteps(obj)
            N=obj.N_samples;S=obj.S_samples;D=obj.D_samples;F=obj.F_samples;M=obj.M_samples;

            I_diff = eye(M) - [zeros(1,M);eye(M-1),zeros(M-1,1)];
            p_M = ones(M,1);
            
            f_fs1 = -obj.Qfsx*obj.vref_x*(obj.ss_time+obj.ds_time)*(I_diff'*p_M);
            f_fs2 = -obj.Qfsy*obj.vref_y*(obj.ss_time+obj.ds_time)*(I_diff'*p_M);
            
            f = [zeros(N,1);f_fs1;zeros(N,1);f_fs2];
        end
        
%         function H = H_obstacle(obj)
%             N=obj.N_samples;S=obj.S_samples;D=obj.D_samples;F=obj.F_samples;M=obj.M_samples;
% 
%             Q_obs_distance = -10;
%             
%             H = Q_obs_distance*[zeros(N,N),zeros(N,M),zeros(N,N),zeros(N,M); ...
%                                 zeros(M,N),  eye(M)  ,zeros(M,N),zeros(M,M); ...
%                                 zeros(N,N),zeros(N,M),zeros(N,N),zeros(N,M); ...
%                                 zeros(M,N),zeros(M,M),zeros(M,N),  eye(M) ];
%         end
%         
%         function f = f_obstacle(obj)
%             N=obj.N_samples;S=obj.S_samples;D=obj.D_samples;F=obj.F_samples;M=obj.M_samples;
%             
%             p_M = ones(M,1);
%             Q_obs_distance = -10;
%             
%             f_fs1 = -Q_obs_distance*obj.obstacle_closest_point(1)*p_M;
%             f_fs2 = -Q_obs_distance*obj.obstacle_closest_point(2)*p_M;
%             
%             f = [zeros(N,1);f_fs1;zeros(N,1);f_fs2];
%         end
        
        function point_in_robot_frame = point_world_to_robot(obj, point_in_world_frame)
            point_in_robot_frame = [cos(obj.theta),-sin(obj.theta);sin(obj.theta),cos(obj.theta)]'*...
                                   (point_in_world_frame - [obj.fxc;obj.fyc]);
        end
        
       
       
       
        
        function cycle(obj,iter)
            N=obj.N_samples;S=obj.S_samples;D=obj.D_samples;F=obj.F_samples;M=obj.M_samples;
            
            % QP options
            options = optimset('Algorithm','interior-point-convex','Display','off');
            
            p = ones(N,1);
            P = tril(ones(N)) * obj.delta;
            
            % Global iteration and iteration since the start of the last single support phase
            obj.iter = iter;
            obj.iterF = mod(obj.iter,F/2);
            
            obj.theta = wrapToPi(obj.theta);
            
            % Current footstep
            obj.footstep = floor(obj.iter/(F/2));

            % Updating time and footsteps
            obj.t_curr = obj.iter*obj.delta;
            rot = [cos(obj.theta),-sin(obj.theta);sin(obj.theta),cos(obj.theta)];
            
            if mod(obj.iter,F/2) == 0 && obj.footstep > 0
                disp(['Changing footstep at ', num2str(obj.t_curr)])
                abs_foot = [obj.fxc;obj.fyc] + rot*[obj.zd(N+1);obj.zd(2*N+M+1)];
                obj.fxc = abs_foot(1);
                obj.fyc = abs_foot(2);
                
                obj.footsteps(:,obj.footstep+1) = [obj.fxc;obj.fyc];
                
                obj.theta = obj.theta + obj.predicted_rotations(2); %pi/16;

                new_rot = [cos(obj.theta),-sin(obj.theta);sin(obj.theta),cos(obj.theta)];
                
                % Change coordinates relative to new support foot
                pos_rel = new_rot'*rot*([obj.x;obj.y] - [obj.zd(N+1);obj.zd(2*N+M+1)]);
                vel_rel = new_rot'*rot*[obj.xd;obj.yd];
                zmp_rel = new_rot'*rot*([obj.zx;obj.zy] - [obj.zd(N+1);obj.zd(2*N+M+1)]);

                obj.x = pos_rel(1);
                obj.y = pos_rel(2);
                obj.xd = vel_rel(1);
                obj.yd = vel_rel(2);
                obj.zx = zmp_rel(1);
                obj.zy = zmp_rel(2);
                
                % If footsteps are preassigned remove past footsteps
                if obj.footsteps_are_preassigned == 1
                    obj.preassigned_footsteps_matrix(1,:) = [];
                end
            end
            
            % Assigning foot rotations
            I_diff_rot = eye(M) - [zeros(1,M);eye(M-1),zeros(M-1,1)];
            p_M = ones(M,1);
            
            if obj.obstacle_is_on
                obj.obstacle_theta = wrapToPi(angdiff(obj.obstacle_theta_abs,obj.theta));
                theta_align = obj.obstacle_theta;
                obj.obstacle_closest_point = obj.point_world_to_robot(obj.obstacle_closest_point_abs);
                obj.obstacle_tangent_point = obj.point_world_to_robot(obj.obstacle_tangent_point_abs);
                d = sqrt((obj.x-obj.obstacle_closest_point(1))^2 + (obj.y-obj.obstacle_closest_point(2))^2);
                theta_gain = wrapToPi(atan2(obj.obstacle_closest_point(2),obj.obstacle_closest_point(1)));
                
                if obj.vref_x < 0
                    if theta_gain > -pi/2 && theta_gain < pi/2
                        ang_gain = 0;
                    else
                        ang_gain = abs(theta_gain) - pi/2;
                    end
                else
                    if theta_gain > -pi/2 && theta_gain < pi/2
                        ang_gain = -abs(theta_gain) + pi/2;
                    else
                        ang_gain = 0;
                    end
                end
            else
                theta_align = 0;
                ang_gain = 0;
                d = 1;
            end
            
            H_rot = obj.Q_omega*(I_diff_rot'*I_diff_rot) + obj.Q_obstacle*eye(M);
            f_rot = -obj.Q_omega*obj.vref_omega*(obj.ss_time+obj.ds_time)*(I_diff_rot'*p_M) - ang_gain*(obj.Q_obstacle*theta_align/d^2)*p_M;
            
            optimal_rotations = quadprog(H_rot,f_rot,[],[],[],[],[],[],[],options);
            
            if obj.footstep == 0
                obj.predicted_rotations = [0,0,optimal_rotations(1:M-1)'];
            else
                obj.predicted_rotations = [0,optimal_rotations(1:M)'];
            end
            
            theta_step_max = pi/16;
            
            for i = 2:M+1
                if wrapToPi(angdiff(obj.predicted_rotations(i), obj.predicted_rotations(i-1))) >= theta_step_max
                    obj.predicted_rotations(i) = wrapToPi(obj.predicted_rotations(i-1) + theta_step_max);
                end
                if wrapToPi(angdiff(obj.predicted_rotations(i), obj.predicted_rotations(i-1))) <= -theta_step_max
                    obj.predicted_rotations(i) = wrapToPi(obj.predicted_rotations(i-1) - theta_step_max);
                end
            end
            
            rot = [cos(obj.theta),-sin(obj.theta);sin(obj.theta),cos(obj.theta)];
            
            % Inequality constraints
            obj.gen_zmp_constraints()
            AQ = obj.A_zmp;
            bQ = obj.b_zmp;
            
            % Add bounding boxes on the footsteps if AFP is on
            if obj.footsteps_are_preassigned == 0
                obj.gen_footstep_constraints()
                AQ = [AQ;obj.A_fs];
                bQ = [bQ;obj.b_fs];
            end
            
            % Add obstacle constraints
            if obj.obstacle_is_on
                obj.gen_obstacle_constraint_on_footsteps()
                AQ = [AQ;obj.A_obstacle];
                bQ = [bQ;obj.b_obstacle];
            end

            % Stability constraint
            %obj.gen_stability_constraint_footsteps()
            obj.gen_stability_constraint()
            
            % Lock the feet in double support
            if obj.iterF > S
                A_lock = zeros(2,2*N+2*M);
                A_lock(1,N+1) = 1;
                A_lock(2,2*N+M+1) = 1;
                b_lock = [obj.zd(N+1);obj.zd(N*2+M+1)];
                obj.Aeq = [obj.Aeq;A_lock];
                obj.beq = [obj.beq;b_lock];
            end
            
            % Lock the first footstep in place
            if obj.footstep == 0 && obj.footsteps_are_preassigned == 0
                A_lock = zeros(2,2*N+2*M);
                A_lock(1,N+1) = 1;
                A_lock(2,2*N+M+1) = 1;
                b_lock = [0;-0.1];
                obj.Aeq = [obj.Aeq;A_lock];
                obj.beq = [obj.beq;b_lock];
            end
            
            % Pre-assigned footsteps
            if obj.footsteps_are_preassigned == 1
                footsteps_remaining = size(obj.preassigned_footsteps_matrix,1);
                predictable_footsteps = min(M,footsteps_remaining);
                
                A_preassigned = zeros(2*predictable_footsteps,2*N+2*M);

                full_A_matrix = eye(M) - [zeros(1,M);eye(M-1),zeros(M-1,1)];
                A_preassigned(1:predictable_footsteps,N+1:N+M) = full_A_matrix(1:predictable_footsteps,1:M);
                A_preassigned(predictable_footsteps+1:2*predictable_footsteps,2*N+M+1:2*(N+M)) = full_A_matrix(1:predictable_footsteps,1:M);
                
                b_preassigned = [obj.preassigned_footsteps_matrix(1:predictable_footsteps,1); ...
                                 obj.preassigned_footsteps_matrix(1:predictable_footsteps,2)];

                obj.Aeq = [obj.Aeq;A_preassigned];
                obj.beq = [obj.beq;b_preassigned];
            end
            
            % Cost function
            obj.gen_cost_function()
            
            % Solver
            obj.zd = quadprog(obj.HQ,obj.fQ,AQ,bQ,obj.Aeq,obj.beq,[],[],[],options);
            
            obj.pred_zx = p * obj.zx + P*obj.zd(1:N);
            obj.pred_zy = p * obj.zy + P*obj.zd(N+M+1:2*N+M);
            obj.pred_fs = [obj.zd(N+1:N+M)';obj.zd(2*N+M+1:2*N+2*M)'];
            
            % Feasibility 
            
            % State update
            ch = cosh(obj.omega*obj.delta);
            sh = sinh(obj.omega*obj.delta);
            A_upd = [ch, sh/obj.omega, 1-ch; obj.omega*sh, ch, -obj.omega*sh; 0, 0, 1];
            B_upd = [obj.delta-sh/obj.omega; 1-ch; obj.delta];
            A_obs = [ch, sh/obj.omega, 1-ch, 0; obj.omega*sh, ch, -obj.omega*sh, obj.delta; 0, 0, 1, 0; 0, 0, 0, 1];
            B_obs = [obj.delta-sh/obj.omega; 1-ch; obj.delta; 0];
            C_obs = [1, 0, 0, 0; 0, 0, 1, 0];
            L = place(A_obs',C_obs', [0.1 0.2 0.12 0.25]);
            
            if obj.iter > 50 && obj.iter < 54
            d_ext_x = 0;
            d_ext_y = 0;
            elseif obj.iter > 74 && obj.iter < 79
            d_ext_x = 0;
            d_ext_y = 0;    
            else
            d_ext_x = 0;
            d_ext_y = 0;
            end
            x_updated = A_upd*[obj.x; obj.xd; obj.zx] + B_upd*obj.zd(1)+[0;obj.delta;0]*d_ext_x;
            y_updated = A_upd*[obj.y; obj.yd; obj.zy] + B_upd*obj.zd(N+M+1)+[0;obj.delta;0]*d_ext_y;
            
            x_hat_updated = A_obs*[obj.x_hat; obj.xd_hat; obj.zx_hat; obj.d_x_hat] + B_obs*obj.zd(1)+  L'*([obj.x;obj.zx]-[obj.x_hat;obj.zx_hat]);
            y_hat_updated = A_obs*[obj.y_hat; obj.yd_hat; obj.zy_hat; obj.d_y_hat] + B_obs*obj.zd(N+M+1)+ L'*([obj.y;obj.zy]-[obj.y_hat;obj.zy_hat]);
  
            
            obj.x = x_updated(1);
            obj.y = y_updated(1);
            obj.xd = x_updated(2);
            obj.yd = y_updated(2);
            obj.zx = x_updated(3);
            obj.zy = y_updated(3);
            
            obj.x_hat = x_hat_updated(1);
            obj.xd_hat = x_hat_updated(2);
            obj.zx_hat = x_hat_updated(3);
            obj.d_x_hat = 0*d_ext_x;
            
            obj.y_hat = y_hat_updated(1);
            obj.yd_hat = y_hat_updated(2);
            obj.zy_hat = y_hat_updated(3);
            obj.d_y_hat = 0*d_ext_y;
            
            
            
            pos_abs = [obj.fxc;obj.fyc] + rot*[obj.x;obj.y];
            vel_abs = rot*[obj.xd;obj.yd];
            zmp_abs = [obj.fxc;obj.fyc] + rot*[obj.zx;obj.zy];
            obj.x_store(iter+1) = pos_abs(1);
            obj.y_store(iter+1) = pos_abs(2);
            obj.xd_store(iter+1) = vel_abs(1);
            obj.yd_store(iter+1) = vel_abs(2);
            obj.zx_store(iter+1) = zmp_abs(1);
            obj.zy_store(iter+1) = zmp_abs(2);
            obj.zxd_store(iter) = obj.zd(1);
            obj.zyd_store(iter) = obj.zd(N+M+1);
                        
            % Clear any obstacle
            obj.A_obstacle = [];
            obj.b_obstacle = [];
        end
        
        function set_plot_limits(obj, xmin, xmax, ymin, ymax)
            obj.plot_limits = [xmin, xmax, ymin, ymax];
        end
        
        function set_plot_options(obj, plot_options)
            obj.plot_options = plot_options;
        end
        
        function plot(obj, iter)
            N=obj.N_samples;S=obj.S_samples;D=obj.D_samples;F=obj.F_samples;M=obj.M_samples;
            
            % Rotate the prediction in the world frame
            pred_to_rotate = [obj.pred_zx , obj.pred_zy]';
            rot = [cos(obj.theta),-sin(obj.theta);sin(obj.theta),cos(obj.theta)];
            pred_rotated = rot*pred_to_rotate;
            
            % Plot the predicted ZMP
            if obj.plot_options.plot_pred_zmp
                plot(obj.fxc+pred_rotated(1,:),obj.fyc+pred_rotated(2,:),'m','LineWidth', 2);
            end
            
            % Plot the CoM
            if obj.plot_options.plot_com
                plot(obj.x_store,obj.y_store,'r','LineWidth', 2);
            end
            
            % Plot the ZMP
            if obj.plot_options.plot_zmp
                plot(obj.zx_store,obj.zy_store,'b','LineWidth', 2);
            end
            
            % Plot the footsteps
%             if obj.plot_options.plot_footsteps
%                 for i = 1 : size(obj.footsteps,2)
%                     plot(obj.footsteps(1,i),obj.footsteps(2,i),'ok');
%                 end
%        
%             end
            
            % Plot the orientation
%             if obj.plot_options.plot_orientation
%                 quiver(obj.x_store(obj.iter),obj.y_store(obj.iter),0.1*cos(obj.theta),0.1*sin(obj.theta),...
%                        'LineWidth',1,'MaxHeadSize',2)
%             end

            % Plot predicted footsteps
            rotated_fs = rot*obj.pred_fs;
            rotated_fs = [[0;0],rotated_fs];
            
            for i = 1 : M+1
%                 Plot predicted footsteps
%                 if obj.plot_options.plot_pred_footsteps
%                     plot(obj.fxc+rotated_fs(1,i),obj.fyc+rotated_fs(2,i),'ok');
%                     hold on
%                 end
                
                % Plot predicted zmp constraints
%                 if obj.plot_options.plot_pred_zmp_constraints
%                     obj.plot_constraint(obj.fxc+rotated_fs(1,i), obj.fyc+rotated_fs(2,i), ...
%                                         obj.w, obj.w, obj.theta+obj.predicted_rotations(i), 'b', 0.1);
%                 
%                 end
                
                % Plot predicted footstep constraints
%                 if obj.plot_options.plot_pred_footstep_constraints && i< M+1
%                     if mod(obj.footstep + i,2) == 0
%                         fs_sign = 1;
%                     else
%                         fs_sign = -1;
%                     end
%                     obj.plot_constraint(obj.fxc+rotated_fs(1,i) - fs_sign*(obj.box_b/2+obj.box_a/2)*sin(obj.theta+obj.predicted_rotations(i)), ...
%                                         obj.fyc+rotated_fs(2,i) + fs_sign*(obj.box_b/2+obj.box_a/2)*cos(obj.theta+obj.predicted_rotations(i)), ...
%                                         obj.box_c*2, obj.box_b-obj.box_a, obj.theta+obj.predicted_rotations(i), [204 255 204]/255, 0.7);
%                 end
            end

            %hold off
            axis equal;
%             xlim([obj.plot_limits(1),obj.plot_limits(2)]); ylim([obj.plot_limits(3),obj.plot_limits(4)]);

            %drawnow
        end
        
        function plot_constraint(obj, xc, yc, xw, yw, theta_c, color, alpha)
            center = [xc xc xc xc; yc yc yc yc];
            vertices = [xw/2 xw/2 -xw/2 -xw/2; +yw/2 -yw/2 -yw/2 +yw/2];
            rotation_matrix = [cos(theta_c), -sin(theta_c); sin(theta_c), cos(theta_c)];
            vertices = rotation_matrix*vertices;
            vertices = center + vertices;
            
            p = patch(vertices(1,:),vertices(2,:),color);
            set(p,'FaceAlpha',alpha,'EdgeColor','k','LineWidth',1,'LineStyle','-');
        end
        
        function make_video(obj, fileName)
            N=obj.N_samples;S=obj.S_samples;D=obj.D_samples;F=obj.F_samples;M=obj.M_samples;
            
            close all;
            
            addpath('functions');

            filePath  = '/home/nickstu/Videos/MPC/';

            fontSizePlot = 14;%
            
            % preallocate
            MFrame(floor(size(obj.x_store,2)/4)) = struct('cdata',[],'colormap',[]);
            
            % Plot
            hfig = figure;
            set(gcf,'Position', [1 1 1024 768]);
            set(gca,'fontsize',fontSizePlot);
            
            samplingTime = obj.delta;
            
            % Linear motion (I)
            xLimits = [-0.1 0.5]; yLimits = [-0.1 0.1];
            
            for i=1:size(obj.x_store,2)
                
                % ZMP
                plot(obj.zx_store(1:i),obj.zy_store(1:i),'b','LineWidth', 2);
                hold on;
                % CoM
                plot(obj.x_store(1:i),obj.y_store(1:i),'r','LineWidth', 2);
                
                grid on;
                hold off;
                
                legend('ZMP', 'CoM', 'Location', 'NorthWest');
                xlabel('x [m]', 'FontSize', fontSizePlot); ylabel('y [m]', 'FontSize', fontSizePlot);
                
                axis equal;
                axis([xLimits yLimits]);
                
                if (mod(i,4)==0)
                    MFrame(i/4) = getframe(hfig);
                end
                
                drawnow
            end
            
            movie2avi(MFrame, [filePath,fileName], 'FPS', 25, 'COLORMAP','jet');
        end
    end
    
end

