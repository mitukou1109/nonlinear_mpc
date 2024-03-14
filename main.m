function main
    f = @(x, u, t) [x(2); (1 - x(1) ^ 2 - x(2) ^ 2) * x(2) - x(1) + u(1)];
    C = @(x, u, t) u(1) ^ 2 + u(2) ^ 2 - 0.5 ^ 2;
    % phi = @(x, t) (x(1) ^ 2 * t + x(2) ^ 2 * t) / 2;
    % L = @(x, u, t) (x(1) ^ 2 + x(2) ^ 2 + u(1) ^ 2) / 2 - 0.01 * u(2);
    % H = @(x, u, lambda, t) L(x, u, t) + lambda.' * f(x, u, t) + u(3) * C(x, u, t);

    DphiDx = @(x, t) [x(1) * t, x(2) * t];
    DHDx = @(x, u, lambda, t) [x(1) + lambda(2) * (-2 * x(1) * x(2) - 1), x(2) + lambda(1) + lambda(2) * (-x(1) ^ 2 - 3 * x(2) ^ 2)];
    DHDu = @(x, u, lambda, t) [u(1) + lambda(2) + 2 * u(1) * u(3), -0.01 + 2 * u(2) * u(3), C(x, u, t)];

    Dt = 0.01;
    T = @(t) 1 - exp(-0.5 * t);
    N = 10;
    zeta = 100;
    GMRESIterations = 2;

    Dtau = @(t) T(t) / N;
    TMax = 20;
    h = 1e-6;

    x = [2; 0];
    u = [0; 0.5; 0];
    t = 0;

    sampleSize = ceil(TMax / Dt);
    tData = zeros(1, sampleSize);
    xData = zeros(numel(x), sampleSize);
    uData = zeros(numel(u), sampleSize);
    FData = zeros(1, sampleSize);

    DUDt = zeros(numel(u), N);

    u = fsolve(@(u) DHDu(x, u, DphiDx(x, t).', t), u);
    U = repmat(u, 1, N);

    for i = 1:sampleSize
        u = U(:, 1);

        Ft = F(U, x, t);

        tData(:, i) = t;
        xData(:, i) = x;
        uData(:, i) = u;
        FData(:, i) = norm(Ft);

        DF = @(W, w, omega) (F(U + h * W, x + h * w, t + h * omega) - Ft) / h;
        g = -zeta * Ft;
        DUDt = fdgmres(DF, g, DUDt, f(x, u, t), 1, GMRESIterations);

        U = U + DUDt * Dt;
        x = x + f(x, u, t) * Dt;
        t = t + Dt;
    end

    %%

    tiledlayout(3, 2)
    nexttile
    plot(tData, xData(1, :))
    % ylim([-1, 2])
    xlabel("Time")
    ylabel("$x_1$", "Interpreter", "latex")
    nexttile
    plot(tData, xData(2, :))
    % ylim([-1, 0.5])
    xlabel("Time")
    ylabel("$x_2$", "Interpreter", "latex")
    nexttile
    plot(tData, uData(1, :))
    % ylim([-0.6, 0.6])
    xlabel("Time")
    ylabel("$u$", "Interpreter", "latex")
    nexttile
    plot(tData, uData(2, :))
    % ylim([0, 0.8])
    xlabel("Time")
    ylabel("$v$", "Interpreter", "latex")
    nexttile
    plot(tData, uData(3, :))
    % ylim([0, 1])
    xlabel("Time")
    ylabel("$\rho$", "Interpreter", "latex")
    nexttile
    plot(tData, FData)
    % ylim([0, 2e-3])
    xlabel("Time")
    ylabel("$\| F \|$", "Interpreter", "latex")

    function [X, LAMBDA] = EulerLagrange(U, x, t)
        X = zeros(numel(x), N);
        LAMBDA = zeros(size(X));

        X(:, 1) = x;

        for j = 1:N - 1
            X(:, j + 1) = X(:, j) + f(X(:, j), U(:, j), t + (j - 1) * Dtau(t)) * Dtau(t);
        end

        LAMBDA(:, N) = DphiDx(X(:, N), t + T(t)).';

        for j = N - 1:-1:1
            LAMBDA(:, j) = LAMBDA(:, j + 1) + DHDx(X(:, j + 1), U(:, j + 1), LAMBDA(:, j + 1), t + (j - 1) * Dtau(t)).' * Dtau(t);
        end

    end

    function F = F(U, x, t)
        F = zeros(size(U));

        [X, LAMBDA] = EulerLagrange(U, x, t);

        for k = 1:size(U, 2)
            F(:, k) = DHDu(X(:, k), U(:, k), LAMBDA(:, k), t + (k - 1) * Dtau(t)).';
        end

        F = reshape(F, [], 1);
    end

end
