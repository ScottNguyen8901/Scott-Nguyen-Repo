function Hh = counted_H(H, C)
% Wraps H so every call increments C.H
    Hh = @(x) local_H(x);
    function y = local_H(x)
        C.H = C.H + 1;
        y = H(x);
    end
end
