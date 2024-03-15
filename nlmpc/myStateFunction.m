function DxDt = myStateFunction(x, u)
    DxDt = [x(2); (1 - x(1) ^ 2 - x(2) ^ 2) * x(2) - x(1) + u(1)];
end
