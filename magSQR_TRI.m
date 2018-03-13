function magSQR_TRI(T_min, T_max, T_step, J, mu, diff_step, bond_dim, log4_N, eps)
    magFun = @(beta) magSQR_TRISub(beta, J, mu, diff_step, bond_dim, log4_N, eps);
    samplesX = T_min:T_step:T_max;
    samplesBeta = 1 ./ samplesX;
    samplesY = arrayfun(magFun, samplesBeta);
    figure(1);
    plot(samplesX, samplesY);
end

function [mag] = magSQR_TRISub(beta, J, mu, diff_step, bond_dim, log4_N, eps)
    logZ1 = partitionSQR_TRI(beta, J, diff_step, mu, bond_dim, log4_N, eps) / beta;
    logZ2 = partitionSQR_TRI(beta, J, 0, mu, bond_dim, log4_N, eps) / beta;
    mag = (logZ1 - logZ2) / diff_step;
end