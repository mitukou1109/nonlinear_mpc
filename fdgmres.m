function W = fdgmres(DF, g, W0, w, omega, iter, tol)

    arguments
        DF
        g
        W0
        w
        omega
        iter
        tol = 1e-6
    end

    r0 = g - DF(W0, w, omega);

    if norm(r0) <= tol
        W = W0;
        return
    end

    V = zeros(numel(W0), iter + 1);
    V(:, 1) = r0 / norm(r0);
    H = zeros(iter + 1, iter);
    x = reshape(W0, [], 1);

    for j = 1:iter
        Avj = DF(reshape(V(:, j), size(W0)), 0, 0);

        for i = 1:j
            H(i, j) = Avj.' * V(:, i);
        end

        s = 0;

        for i = 1:j
            s = s + H(i, j) * V(:, i);
        end

        vHat = Avj - s;

        H(j + 1, j) = norm(vHat);
        V(:, j + 1) = vHat / H(j + 1, j);
    end

    y = lsqminnorm(H, norm(r0) * eye(iter + 1, 1));
    x = x + V(:, 1:iter) * y;
    W = reshape(x, size(W0));
end
