classdef rigidity_dlmpc_controller < Controller
    properties
        N
        model %Model = AUV(ones(11,1),3,zeros(6,1),zeros(3,1))
        auxiliary_controller 
        weights
        upperbound
        lowerbound
        Adj
    end
    
    methods
        function obj = rigidity_dlmpc_controller(N,model,controller,weights,upperbound,lowerbound,Adj)
            obj.N = N;
            obj.model = model;
            obj.auxiliary_controller = controller;
            obj.weights = weights;
            obj.calc_upperbound(upperbound);
            obj.calc_lowerbound(lowerbound);
            obj.Adj = Adj;
        end
        
        % 上界
        function obj = calc_upperbound(obj,upperbound)
            if ~isempty(upperbound)
                m = length(upperbound);
                U_max = zeros(m*obj.N,1);
                for i=1:1:obj.N
                    U_max(1+m*(i-1):m*i,1) = upperbound;
                end
                obj.upperbound = U_max;
            else
                obj.upperbound = [];
            end
        end
        
        % 下界
        function obj = calc_lowerbound(obj,lowerbound)
            if ~isempty(lowerbound)
                m = length(lowerbound);
                U_min = zeros(m*obj.N,1);
                for i=1:1:obj.N
                    U_min(1+m*(i-1):m*i,1) = lowerbound;
                end
                obj.lowerbound = U_min;
            else
                obj.lowerbound = [];
            end
        end
        
        % 代价函数cost function
        function cost = dlmpc_cost(obj,u,ref,id,traj,state,state_pred,dt)
            % 被控机器人初始状态
            X0 = state(:,id);
            % 预测区间
            Hp = obj.N;
            % 权重
            L = obj.weights{1};
            Q = obj.weights{2};
            R = obj.weights{3};
            Qf = obj.weights{4};
            % 输入u的长度
            nu = length(u);
            nu = nu/Hp;
            % 状态x的长度
            nx = length(X0);
            Hu = Hp;
            % 初始化U X
            U = zeros(nu,Hp);
            X = zeros(nx,Hp);
           
            psi_r = zeros(Hp,1);
            ur = zeros(Hp,1);
            vr = zeros(Hp,1);
            rr = zeros(Hp,1);
         
            for i=1:1:Hp
                psi_r(i,1) = traj(3,i);
                ur(i,1) = traj(4,i);
                vr(i,1) = traj(5,i);
                rr(i,1) = traj(6,i);
            end

            % 把列向量u划分为nu行Hp列的矩阵，每一列代表每一步预测
            % partition of u
            for i=1:1:Hp
                for j=1:1:nu
                    U(j,i) = u((i-1)*nu+j,1);
                end
            end
            
            % AUV运动学&动力学计算每次预测的输出所对应的状态
            Xplus = obj.model.dynamics_discrete( X0,U(:,1),dt );
            X(:,1) = Xplus;
            for i=2:1:Hp
                Xplus = obj.model.dynamics_discrete( Xplus,U(:,i),dt );
                X(:,i)= Xplus;
            end
            
            % 计算η~
            q_tilde_norm = zeros(length(obj.Adj),Hp);
            % 对预测区间内不同时刻进行计算
            for i=1:Hp
                eta = X(1:3,i);
                q = eta(1:2);
                % 对临近的机器人进行计算
                for j = 1:length(obj.Adj)
                    if(obj.Adj(id,j))
                        eta_j = state_pred(1:3,i,j);
                        q_j = eta_j(1:2);
                        q_tilde = q_j - q;
                        q_tilde_norm(j,i) = norm(q_tilde);
                    end
                end
            end
            % q~与期望值的误差
            error = q_tilde_norm - ref(id,:)';
            cost = 0;
         
            for i=1:1:Hp
                Xr = [ur(i,1);vr(i,1);rr(i,1)];
                v = X(4:6,i);
                psi = X(3,i);
                eta_dot = Rot(psi)*v;
                cost = cost + (psi_r(i)-psi)'*L(1)*(psi_r(i)-psi);
                cost = cost + (eta_dot-Xr)'*L(2:4,2:4)*(eta_dot-Xr);
                if(i<Hp)
                    cost = cost + error(:,i)'*Q*error(:,i); 
                end
            end
           
            % 输出u
            for i=1:1:Hu
                cost = cost + U(:,i)'*R*U(:,i);
            end
            cost = cost + error(:,Hp)'*Qf*error(:,Hp);
        end

        % 最优化算法的非线性约束
        % 使MPC的V函数小于反步控制的V函数
        function [y,ceq] = dlmpc_contraints(obj,u,ref,traj,id,state)
            ceq=[];
            Hp = obj.N;
            nu = length(u); %[U0,U1,..U_N-1]
            nu = nu/Hp; % Ui=[Fu_i, Fv_i, Fr_i]'
            U = zeros(nu,Hp);
            
            % 实际机器人的状态η,v，φ
            eta = state(1:3,id); % [x;y;psi]
            vel = state(4:6,id); % [u;v;r]
            psi = state(3,id);
            q   = eta(1:2); 
            % 旋转矩阵R(φ)
            Rpsi = Rot(psi);

            % 求η'
            eta_dot=Rpsi*vel; % [x';y';psi']
            q_dot = eta_dot(1:2);

            Ka = obj.auxiliary_controller.Ka;
            Kv = obj.auxiliary_controller.Kv;
            K_psi_v = obj.auxiliary_controller.K_psi_v;
            K_psi_a = obj.auxiliary_controller.K_psi_a;
            
            M = diag([obj.model.DeducedCoef]);
            M_star = Rpsi*M*Rpsi';

            % Coriolis and centripetal matrix C(v) 

            C = [0            ,-M(2,2)*vel(3),0;
                 M(1,1)*vel(3),             0,0;
                 0            ,             0,0];
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
                    q_tilde(:,i) = q - q_j; % η~
                    z(i,:) =  q_tilde(1,i)^2 + q_tilde(2,i)^2 - ref(id,i)^2; % z_ij
                    q_tilde_dot(:,i) = q_dot - q_dot_j;% η~'
                end
            end

            I2 = ones(2,2);
            for i = 1:length(z)
                temp(:,i) = ( z(i)*I2 + 2*q_tilde(:,i)*q_tilde(:,i)' )*q_tilde_dot(:,i);
            end

            RTz = q_tilde*z;
            v_d = traj(4:5,1);
            v_d_dot = traj(7:8,1);
            q_f_dot = -Kv*RTz + v_d;
            S = q_dot - q_f_dot;
            q_f_dot_dot = -Kv*sum(temp,2) + v_d_dot;
            
            Rtau_bs = (-Ka*S + C_star(1:2,1:2)*q_f_dot + D_star(1:2,1:2)*q_dot + M_star(1:2,1:2)*q_f_dot_dot - RTz);

            Vdot_bs(1,1) = RTz'*q_f_dot + S'*(Rtau_bs - C_star(1:2,1:2)*q_f_dot - D_star(1:2,1:2)*q_dot - M_star(1:2,1:2)*q_f_dot_dot + RTz);
            
            psi_r = traj(3,1);
            r_r = traj(6,1);
            r_r_dot = 0;
            Mx = obj.model.DeducedCoef(1);
            My = obj.model.DeducedCoef(2);
            Mpsi = obj.model.DeducedCoef(3);

            A_r = ((Mx - My)*vel(1)*vel(2) - D(3,3)*vel(3))/Mpsi;
            r_f = r_r + (psi_r - psi)*K_psi_v;
            r_f_dot = r_r_dot + (r_r - vel(3))*K_psi_v;
            s_r = vel(3) - r_f;
            tau_bs(3) = (-A_r + r_f_dot + psi_r - psi - K_psi_a*s_r)*Mpsi;
            Vdot_bs(2,1) = s_r*(tau_bs(3)/Mpsi + A_r - r_f_dot - psi_r + psi) +(psi_r - psi)*(r_r - r_f);
            
            for i=1:1:Hp
                for j=1:1:nu
                    U(j,i)=u((i-1)*nu+j,1);
                end
            end

            RTau=Rpsi*U(:,1);
            RTau(3)=[];
            Vdot(1,1) = RTz'*q_f_dot + S'*(RTau - C_star(1:2,1:2)*q_f_dot - D_star(1:2,1:2)*q_dot - M_star(1:2,1:2)*q_f_dot_dot + RTz);
            Tau=U(3,1);
            Vdot(2,1) = s_r*(Tau/Mpsi + A_r - r_f_dot - psi_r + psi) +(psi_r - psi)*(r_r - r_f);
            y = Vdot - Vdot_bs;
        end
        
        % general control
        function [u, X] = calc_control(obj,ref,id,traj,state,state_pred,u0,dt)
            options = optimset('Algorithm','sqp');
            % ,'Display','none');
            % Calculate the MPC control input 
            u = fmincon(@(u) obj.dlmpc_cost(u,ref,id,traj,state,state_pred,dt),u0,[],[],[],[],...
                obj.lowerbound, obj.upperbound, ...
                @(u) obj.dlmpc_contraints(u,ref,traj,id,state),options);
            % u = fmincon(@(u) obj.dlmpc_cost(u,ref,id,traj,state,state_pred,dt),u0,[],[],[],[],...
            %     obj.lowerbound, obj.upperbound, ...
            %     [],options);
            % Construct the broadcast information
            Hp = obj.N;
            Xplus = obj.model.dynamics_discrete(state(:,id), u(1:3),dt);
            X(:,1) = Xplus;

            for i=2:1:Hp
                Xplus = obj.model.dynamics_discrete(Xplus,u(3*i-2:3*i),dt);
                X(:,i)= Xplus;
            end
        end

        % 以反步控制的输入量为初始值进行最优化计算
        function u0 = calc_initial_guess(obj,ref,traj,id,state,dt)
            Hp = obj.N;
            nu = length(obj.model.U);
            u0 = zeros(nu*Hp,1);
            state_ = state;
            for j = 1:1:Hp
                U_auxiliary = obj.auxiliary_controller.calc_control(ref,traj,id,state_);
                u0(nu*(j-1)+1:nu*j) = U_auxiliary;
                state_(:,id) = obj.model.dynamics_discrete(state_(:,id),U_auxiliary,dt);
            end
        end
        
        function get_params(obj)
            fprintf('weights: \n')
            disp(obj.weights)
        end
    end
end


