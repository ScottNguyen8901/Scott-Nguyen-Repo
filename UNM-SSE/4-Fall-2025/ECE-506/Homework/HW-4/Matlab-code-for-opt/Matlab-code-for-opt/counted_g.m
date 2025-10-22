function gh = counted_g(g, C)
% Wraps g so every call increments C.g
    gh = @(x) local_g(x);
    function y = local_g(x)
        C.g = C.g + 1;
        y = g(x);
    end
end
