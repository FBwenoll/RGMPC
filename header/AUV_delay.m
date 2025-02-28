classdef AUV_delay < handle %Model
    properties
        Coef = [116.0; 13.1;-167.6;-477.2;-15.9;26.9;35.8;3.5;241.3;503.8;76.9]
        Ndof = 3;
        DeducedCoef = [283.6; 593.2; 29.0];
        X
        U
        U_delayed
        alpha
    end
    methods
        function obj = AUV_delay(coef,ndof,X0,U0,dt)
            if nargin ==  5
                obj.Coef = coef;
                obj.Ndof = ndof;
                obj = calc_deduced_coef(obj);
                obj.X = X0;
                obj.U = U0;
                obj.U_delayed = U0; % 初始延迟控制输入等于初始控制输入

                obj.alpha = dt / (0.5 + dt); % 计算 alpha
            else
                errID = 'myComponent:inputError';
                msg = 'Five input arguments (coef, ndof, X0, U0, dt) should be provided.';
                baseException = MException(errID,msg);
                throw(baseException);
            end
        end
        
        function obj = calc_deduced_coef(obj)
            if obj.Ndof == 3
                assert(length(obj.Coef) == 11,'Number of model coeffients is incorrect: For 3 DOF AUV model we need 11 parameters.')
                m = obj.Coef(1);
                Iz = obj.Coef(2);
                X_udot = obj.Coef(3);
                Y_vdot = obj.Coef(4);
                N_rdot = obj.Coef(5);
                Mx = m-X_udot;
                My = m-Y_vdot;
                Mpsi = Iz-N_rdot;
                obj.DeducedCoef = [Mx;My;Mpsi];
            else
                errID = 'myComponent:inputError';
                msg = 'The Ndof is currently not supported';
                baseException = MException(errID,msg);
                throw(baseException);
            end
        end
        
        function Xplus = dynamics_discrete(obj,X,U,dt)
            Mx = obj.DeducedCoef(1);
            My = obj.DeducedCoef(2);
            Mpsi = obj.DeducedCoef(3);
            Xu = obj.Coef(6);
            Yv = obj.Coef(7);
            Nr = obj.Coef(8);
            Du = obj.Coef(9);
            Dv = obj.Coef(10);
            Dr = obj.Coef(11);

            x = X(1);
            y = X(2);
            psi = X(3);
            u = X(4);
            v = X(5);
            r = X(6);

            Fu = U(1);
            Fv = U(2);
            Fr = U(3);

            x_dot = u * cos(psi) - v * sin(psi);
            y_dot = u * sin(psi) + v * cos(psi);
            psi_dot = r;
            u_dot = (My / Mx) * v * r - (Xu / Mx) * u - (Du / Mx) * u * abs(u) + Fu / Mx;
            v_dot = -(Mx / My) * u * r - (Yv / My) * v - (Dv / My) * v * abs(v) + Fv / My;
            r_dot = ((Mx - My) / Mpsi) * u * v - (Nr / Mpsi) * r - (Dr / Mpsi) * r * abs(r) + Fr / Mpsi;

            X_dot = [x_dot; y_dot; psi_dot; u_dot; v_dot; r_dot];

            Xplus = X + X_dot * dt;
        end
        
        function obj = advance(obj,U,W,dt)
            if obj.Ndof == 3
                % 计算延迟控制输入
                obj.U_delayed = obj.alpha * U + (1 - obj.alpha) * obj.U_delayed;
                
                disturbed_U = obj.U_delayed + W;
                obj.X = obj.dynamics_discrete(obj.X,disturbed_U,dt);
                obj.U = U;
            else
                errID = 'myComponent:inputError';
                msg = 'The Ndof is currently not supported';
                baseException = MException(errID,msg);
                throw(baseException);
            end
        end

        % 其他方法保持不变...
    end
end
