function [A, Bmv] = myStateJacobian(x, ~)
    A = [0, 1; -2 * x(1) * x(2) - 1, 1 - x(1) ^ 2 - 3 * x(2)];
    Bmv = [0, 0; 1, 0];
end
