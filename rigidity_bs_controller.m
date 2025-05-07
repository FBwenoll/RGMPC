classdef rigidity_bs_controller < Auxuliary_Controller
    properties
        Kv
        Ka
        K_psi_v
        K_psi_a
        coef
        Adj
        upperbound
        lowerbound
    end
    
    methods
        function obj = rigidity_bs_controller(k_bs,coef,upperbound,lowerbound,Adj)
            obj.Kv = k_bs.v;
            obj.Ka = k_bs.a;
            obj.K_psi_v = k_bs.psi_v;
            obj.K_psi_a = k_bs.psi_a;
            obj.coef = coef;
            obj.Adj = Adj;
            obj.upperbound = upperbound;
            obj.lowerbound = lowerbound;
        end
        
        function [U,Vdot] = calc_control(obj,ref,traj,id,state)
            m = obj.coef(1);
            Iz = obj.coef(2);
            X_udot = obj.coef(3);
            Y_vdot = obj.coef(4);
            N_rdot =  obj.coef(5);
            
            Mx=m-X_udot;
            My=m-Y_vdot;
            Mpsi=Iz-N_rdot;
           
            % 实际机器人的状态η,v，φ
            eta = state(1:3,id); % [x;y;psi]
            vel = state(4:6,id); % [u;v;r]
            psi = state(3,id);
            q   = eta(1:2);
            % 旋转矩阵R(φ)
            Rpsi = Rot(psi);

            % 求η'
            eta_dot=Rpsi*vel; % [x';y';psi']
            q_dot  = eta_dot(1:2);
            % mass matrix
            M = diag([Mx;My;Mpsi]);
            M_star = Rpsi*M*Rpsi';

            % Coriolis and centripetal matrix C(v) 
            C = [0        ,-My*vel(3),0;
                 Mx*vel(3),         0,0;
                 0        ,         0,0];
            psi_dot=eta_dot(3); % φ'
            S_psi_dot = [0      ,-psi_dot,0;
                         psi_dot,0       ,0;
                         0      ,0       ,0];
            C_star = Rpsi*(C-M*S_psi_dot)*Rpsi';

            % damping matrix D(v)
            D = DV(vel);
            D_star = Rpsi*D*Rpsi';

            % 计算参数
            q_tilde = zeros(2,length(obj.Adj));
            z = zeros(length(q_tilde),1);
            q_tilde_dot = zeros(2,length(obj.Adj));
            for i = 1:length(obj.Adj)
                if(obj.Adj(id,i))
                    % 计算控制量中的各参数
                    eta_j = state(1:3,i); % η_j
                    vel_j = state(4:6,i); % v_j
                    Rpsi_j = Rot(eta_j(3));    % R_j
                    eta_dot_j = Rpsi_j*vel_j; % η'_j
                    q_dot_j = eta_dot_j(1:2);
                    q_j = eta_j(1:2);
                    q_tilde(:,i) = q - q_j; % q~
                    z(i,:) =  q_tilde(1,i)^2 + q_tilde(2,i)^2 - ref(id,i)^2; % z_ij
                    q_tilde_dot(:,i) = q_dot - q_dot_j;
                end
            end

            %   Backstepping
            I2 = ones(2,2);
            for i = 1:length(z)
                temp(:,i) = ( z(i)*I2 + 2*q_tilde(:,i)*q_tilde(:,i)' )*q_tilde_dot(:,i);
            end

            RTz = q_tilde*z;
                 
            q_f_dot = -obj.Kv*RTz+traj(2:3)';  
            S = q_dot - q_f_dot;
            
            q_f_dot_dot = -obj.Kv*sum(temp,2)+traj(5:6)';
            
            RU = -obj.Ka*S + C_star(1:2,1:2)*q_f_dot + D_star(1:2)*q_dot + M_star(1:2,1:2)*q_f_dot_dot - RTz;
            Vdot(1) = RTz'*q_f_dot + S'*(RU - C_star(1:2,1:2)*q_f_dot - D_star(1:2,1:2)*q_dot - M_star(1:2,1:2)*q_f_dot_dot + RTz);
            U = Rpsi'*[RU(1);RU(2);0];
            U = [U(1);U(2)];
            
            
            % psi闭环
            psi_r = traj(1);
            r_r = traj(4);
            r_r_dot = 0;
            u = vel(1);
            v = vel(2);
            r = vel(3);
            A_r = ((Mx - My)*u*v - D(3,3)*r)/Mpsi;
            r_f = r_r + (psi_r - psi)*obj.K_psi_v;
            r_f_dot = r_r_dot + (r_r - r)*obj.K_psi_v;
            s_r = r - r_f;
            U(3) = (-A_r + r_f_dot + psi_r - psi - obj.K_psi_a*s_r)*Mpsi;
            Vdot(2) = s_r*(U(3)/Mpsi + A_r - r_f_dot - psi_r + psi) +(psi_r - psi)*(r_r - r_f);
            
            U(U>obj.upperbound) = obj.upperbound;
            U(U<obj.lowerbound) = obj.lowerbound;
        end
            
        function get_params(obj)
            fprintf('K: \n');
        end
        
        function Vdot = lyapunov_derivative(obj,ref,psi_traj,id,state)
            [~,Vdot] = calc_control(obj,ref,psi_traj,id,state);          
        end
    end
end

