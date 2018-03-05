function freeEnergySQR(beta_min, beta_max, step_size, bond_dim, log4_N)

    partitionFun = @(beta) partitionSQR(beta, bond_dim, log4_N) / beta;
    samplesX = beta_min:step_size:beta_max;
    samplesY = arrayfun(partitionFun, samplesX);
    figure(1);
    plot(samplesX, samplesY);
end