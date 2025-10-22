function fh = counted_f(f, C)
% Wraps f so every call increments C.f
    fh = @(x) local_f(x);
    function y = local_f(x)
        C.f = C.f + 1;
        y = f(x);
    end
end
