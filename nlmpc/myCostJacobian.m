function [G, Gmv, Ge] = myCostJacobian(X, U, ~, data)
    G = X(2:end, :);
    Gmv = [U(1:end - 1, 1), -0.01 * ones(data.PredictionHorizon, 1)];
    Ge = 0;
end
