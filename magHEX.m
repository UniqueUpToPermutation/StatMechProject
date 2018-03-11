function magHEX(T_min, T_max, T_step, J, diff_step, bond_dim, log4_N, eps)
    magFun = @(beta) magHEXSub(beta, J, diff_step, bond_dim, log4_N, eps);
    samplesX = T_min:T_step:T_max;
    samplesBeta = 1 ./ samplesX;
    samplesY = arrayfun(magFun, samplesBeta);
    figure(1);
    plot(samplesX, samplesY);
end

function [mag] = magHEXSub(beta, J, diff_step, bond_dim, log4_N, eps)
    logZ1 = partitionHEX(beta, J, diff_step, bond_dim, log4_N, eps) / beta;
    logZ2 = partitionHEX(beta, J, 0, bond_dim, log4_N, eps) / beta;
    mag = (logZ1 - logZ2) / diff_step;
end