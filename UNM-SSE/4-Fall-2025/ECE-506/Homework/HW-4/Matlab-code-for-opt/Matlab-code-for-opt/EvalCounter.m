classdef EvalCounter < handle
    properties
        f = 0;   % # of function evaluations
        g = 0;   % # of gradient evaluations
        H = 0;   % # of Hessian evaluations
    end
    methods
        function reset(self)
            self.f = 0; self.g = 0; self.H = 0;
        end
    end
end
