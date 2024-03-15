function J = myCostFunction(X, U, ~, data)
    J = 0;

    for i = 2:data.PredictionHorizon
        J = J + (X(i, 1) ^ 2 + X(i, 2) ^ 2 + U(i - 1, 1) ^ 2) / 2 - 0.01 * U(i - 1, 2);
    end

    J = J + (X(end, 1) ^ 2 + X(end, 2) ^ 2) / 2;

end
