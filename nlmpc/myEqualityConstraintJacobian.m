function [Geq, Gmv] = myEqualityConstraintJacobian(~, U, data)
    U = U(1:data.PredictionHorizon, :);
    Geq = zeros(data.PredictionHorizon, 2, data.PredictionHorizon);
    Gmv = permute(cat(3, diag(2 * U(:, 1)), diag(2 * U(:, 2))), [1, 3, 2]);
end
