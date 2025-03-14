classdef linear2D < Trajectory2D
    properties
        slope
        intercept
        x_initial
    end
    methods
        function obj = linear2D(m, b, x0, t0)
            if nargin == 4
                obj.slope = m;
                obj.intercept = b;
                obj.x_initial = x0;
                obj.t = t0;
                obj = init(obj, t0);
            else
                errID = 'myComponent:inputError';
                msg = 'Four input arguments (slope, intercept, x_initial, t0) should be provided.';
                baseException = MException(errID, msg);
                throw(baseException);
            end
        end
        
        function xR = Xr(obj, t)
            xR = obj.x_initial + t; % 调整 x 坐标的变化速率
        end
        
        function xRd = Xrdot(obj, t)
            xRd = 0.05; % dx/dt 的变化速率调整为 0.1
        end
        
        function xRdd = Xrddot(obj, t)
            xRdd = 0; % d^2x/dt^2 仍为零
        end
        
        function xRddd = Xrdddot(obj, t)
            xRddd = 0; % d^3x/dt^3 仍为零
        end
        
        function yR = Yr(obj, t)
            x = obj.Xr(t);
            yR = obj.slope * x + obj.intercept;
        end
        
        function yRd = Yrdot(obj, t)
            yRd = obj.slope * obj.Xrdot(t); % dy/dt = slope * dx/dt
        end
        
        function yRdd = Yrddot(obj, t)
            yRdd = 0; % d^2y/dt^2 仍为零
        end
        
        function yRddd = Yrdddot(obj, t)
            yRddd = 0; % d^3y/dt^3 仍为零
        end
    end
end
