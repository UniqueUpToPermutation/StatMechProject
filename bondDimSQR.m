function bondDimSQR(T_min, T_max, T_step, J, bond_dim, log4_N, eps)
    magFun = @(beta) magSQRSub(beta, J, bond_dim, log4_N, eps);
    samplesX = T_min:T_step:T_max;
    samplesBeta = 1 ./ samplesX;
    samplesBetaCell = num2cell(samplesBeta');
    samplesY = cellfun(magFun, samplesBetaCell, 'un', 0);
    samples = cell2mat(samplesY)';
    figure(1);
    imagesc(samplesX, 1:size(samples, 1), samples);
    colorbar;
    set(gca,'YDir','normal')
end

function [bond_dims] = magSQRSub(beta, J, bond_dim, log4_N, eps)
    [~, aux] = partitionSQR(beta, J, 0, bond_dim, log4_N, eps);
    bond_dims = aux.bond_dims;
end