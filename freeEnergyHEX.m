function freeEnergyHEX(T_min, T_max, T_step, J, h, bond_dim, log4_N, eps)
    partitionFun = @(beta) partitionSQR(beta, J, h, bond_dim, log4_N, eps) / beta;
    samplesX = T_min:T_step:T_max;
    samplesBeta = 1 ./ samplesX;
    samplesY = arrayfun(partitionFun, samplesBeta);
    figure(1);
    plot(samplesX, samplesY);
end