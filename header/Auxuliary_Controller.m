classdef (Abstract) Auxuliary_Controller < Controller
    methods(Abstract)
         Vdot = lyapunov_derivative(ref,id,state)
    end 
end