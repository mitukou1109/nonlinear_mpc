function ceq = myEqualityConstraintFunction(~, U, data)
    U = U(1:data.PredictionHorizon, :);
    ceq = U(:, 1) .^ 2 + U(:, 2) .^ 2 - 0.5 ^ 2;
end
