function W = fdgmres(DF, g, W0, w, omega, iter)

    arguments
        DF
        g
        W0
        w
        omega
        iter
    end

    W = W0;
    x = reshape(W0, [], 1);
    m = numel(W0);

    for k = 1:iter
        V = zeros(m, m + 1);
        H = zeros(m + 1, m);

        r0 = g - DF(W, w, omega);
        V(:, 1) = r0 / norm(r0);

        for j = 1:m
            Avj = DF(reshape(V(:, j), size(W0)), 0, 0);

            for i = 1:j
                H(i, j) = Avj.' * V(:, i);
            end

            s = zeros(size(Avj));

            for i = 1:j
                s = s + H(i, j) * V(:, i);
            end

            vHat = Avj - s;

            H(j + 1, j) = norm(vHat);
            V(:, j + 1) = vHat / H(j + 1, j);
        end

        y = lsqminnorm(H, norm(r0) * eye(m + 1, 1));
        x = x + V(:, 1:m) * y;
        W = reshape(x, size(W0));
    end

end
