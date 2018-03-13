function magSQR(T_min, T_max, T_step, J, diff_step, bond_dim, log4_N, eps)
    magFun = @(beta) magSQRSub(beta, J, diff_step, bond_dim, log4_N, eps);
    samplesX = T_min:T_step:T_max;
    samplesBeta = 1 ./ samplesX;
    samplesY = arrayfun(magFun, samplesBeta);
    figure(1);
    plot(samplesX, samplesY);
end

function [mag] = magSQRSub(beta, J, diff_step, bond_dim, log4_N, eps)
    [logZ1, ~] = partitionSQR(beta, J, diff_step, bond_dim, log4_N, eps);
    f1 = logZ1 / beta;
    [logZ2, ~] = partitionSQR(beta, J, 0, bond_dim, log4_N, eps);
    f2 = logZ2 / beta;
    mag = (f1 - f2) / diff_step;
end