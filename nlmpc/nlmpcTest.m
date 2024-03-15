Dt = 0.01;
T = @(t) 1 - exp(-0.5 * t);
N = 10;

TMax = 20;
t = 0;
x = [2; 0];
u = [0; 0.5];

nlobj = nlmpc(numel(x), numel(x), numel(u));
nlobj.PredictionHorizon = N;
nlobj.ControlHorizon = nlobj.PredictionHorizon;
nlobj.Model.StateFcn = @myStateFunction;
nlobj.Optimization.CustomCostFcn = @myCostFunction;
nlobj.Optimization.CustomEqConFcn = @myEqualityConstraintFunction;
nlobj.Jacobian.StateFcn = @myStateJacobian;
nlobj.Jacobian.CustomCostFcn = @myCostJacobian;
nlobj.Jacobian.CustomEqConFcn = @myEqualityConstraintJacobian;

validateFcns(nlobj, x, u);

sampleSize = ceil(TMax / Dt);
tData = zeros(1, sampleSize);
xData = zeros(numel(x), sampleSize);
uData = zeros(numel(u), sampleSize);

for i = 1:sampleSize
    tData(:, i) = t;
    xData(:, i) = x;
    uData(:, i) = u;

    if T(t) / N > 0
        nlobj.Ts = T(t) / N;
        u = nlmpcmove(nlobj, x, u);
    end

    x = x + myStateFunction(x, u) * Dt;
    t = t + Dt;
end

%%

tiledlayout(2, 2)
nexttile
plot(tData, xData(1, :))
ylim([-1, 2])
xlabel("Time")
ylabel("$x_1$", "Interpreter", "latex")
nexttile
plot(tData, xData(2, :))
ylim([-1, 0.5])
xlabel("Time")
ylabel("$x_2$", "Interpreter", "latex")
nexttile
plot(tData, uData(1, :))
ylim([-0.6, 0.6])
xlabel("Time")
ylabel("$u$", "Interpreter", "latex")
nexttile
plot(tData, uData(2, :))
ylim([0, 0.8])
xlabel("Time")
ylabel("$v$", "Interpreter", "latex")
