function varSQR_TRI_mu(K, mu_start, mu_stop, mu_step, diff_step, bond_dim, log4_N, eps)
    % magFun = @(mu) magSQR_TRISub(K, mu, diff_step, bond_dim, log4_N, eps);
    secFun = @(mu) secSQR_TRISub(K, mu, diff_step, bond_dim, log4_N, eps);
    samplesX = mu_start:mu_step:mu_stop;
    % samplesMag = arrayfun(magFun, samplesX);
    samplesSecondMom = arrayfun(secFun, samplesX);
    % samplesY = samplesSecondMom - samplesMag .^ 2;
    samplesY = samplesSecondMom;
    figure(1);
    plot(samplesX, samplesY);
end

%{
function [mag] = magSQR_TRISub(K, mu, diff_step, bond_dim, log4_N, eps)
    logZ1 = partitionSQR_TRI(1, K, diff_step, mu, bond_dim, log4_N, eps) / 1;
    logZ2 = partitionSQR_TRI(1, K, 0, mu, bond_dim, log4_N, eps) / 1;
    mag = (logZ1 - logZ2) / diff_step;
end
%}

function [mag] = secSQR_TRISub(K, mu, diff_step, bond_dim, log4_N, eps)
    logZ1 = partitionSQR_TRI(1, K, 0, mu + diff_step, bond_dim, log4_N, eps) / 1;
    logZ2 = partitionSQR_TRI(1, K, 0, mu, bond_dim, log4_N, eps) / 1;
    mag = (logZ1 - logZ2) / diff_step;
end