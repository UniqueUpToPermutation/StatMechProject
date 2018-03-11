% Compute the log partition function for a honeycomb lattice
function [logZ_perSite] = partitionHEX(beta, J, h, bond_dim, log4_N, eps)
    T = tensorHEX(beta, J, h);
    
    if nargin < 4
        eps = 1E-2;
    end
    
    fprintf('\nBeta = %f\n', beta);
    fprintf('-----------------------------------------\n');
    
    logZ_perSite = 0;
    log2_N = 2 * log4_N;
    factor = 1;
   
    for log2_n=log2_N:-1:2 % Log4 number of lattice sites
        M = tensorHEXJoin(T);
        % Split the tensor T into two copies of S
        [S, sigma1] = tensorHEXSplit(M, bond_dim, eps);
        % Contract four copies of S around a loop to form the new T
        T = loopContractHEX(S);
        
        % Normalize by sigma1
        factor = factor / 3;
        T = T / sigma1;
        logZ_perSite = logZ_perSite + factor * log(sigma1);
    end
    
    fprintf('-----------------------------------------\n');
end

% Perform a loop contraction of three S tensors
function [T] = loopContractHEX(S)
    % T'_{ijk} = sum_{i'j'k'} S_{i'j'i} S_{j'k'j} S_{k'i'k}
    T1 = ttt(S, S, 1, 1);
    T = ttt(T1, S, [1, 3], [2, 1]);
end

function [delta] = accuracyCheck(S, T)
    % T_{ijkl} = sum_m S_{ijm} S_{klm}
    Tapprox = ttt(S, S, 3, 3);
    diff = T - Tapprox;
    delta = norm(diff);
    nrm = norm(T);
    delta = delta / nrm;
end

% Perform a join of the T tensor
function [M] = tensorHEXJoin(T)
    % M_{ijkl} = sum_m T_{ijm} T_{klm}
    M = ttt(T, T, 1, 1);
end

% Perform a split of the M tensor
function [S, sigma1] = tensorHEXSplit(M, bond_dim, eps)
    i_dim = size(M, 1);
    j_dim = size(M, 2);
    k_dim = size(M, 3);
    l_dim = size(M, 4);
    assert(i_dim == j_dim && j_dim == k_dim && k_dim == l_dim);
    
    % Turn the tensor T into a matrix for splitting
    T_ = tenmat(M, [1, 2], [3, 4]);
    
    % Take SVD of T
    [U, sigmas, ~] = svd(T_.data);
    
    diags_ = diag(sigmas);
    diags_ = diags_(diags_ >= eps * diags_(1));
    rank_ = length(diags_);
    
    % Truncate using the bond_dimension
    bond_dim = min(bond_dim, rank_);
    U_ = U(:,1:bond_dim);
    diags_ = sigmas(1:bond_dim, 1:bond_dim);
    S_ = U_ * sqrt(diags_);
    
    % Reshape into split tensor S
    S = tensor(S_);
    S = reshape(S, [i_dim, j_dim, bond_dim]);
    
    % Check that the approximation is accurate
    err = accuracyCheck(S, M);
    fprintf('Relative Approximation Error: %f\t[Bond Dimension = %i]\n', ...
        err, bond_dim);
    
    sigma1 = diags_(1, 1);
end

% Setup initial T-tensor
function [T] = tensorHEX(beta, J, h)
    d = 3;

    bond_matrix_J = [exp(J * beta), exp(-J * beta); exp(-J * beta), exp(J * beta)];
    bond_matrix_h = [exp(2 / d * h * beta), 1; 1, exp(-2 / d * h * beta)];
    bond_matrix = bond_matrix_J .* bond_matrix_h;
    [U, Sigma, V] = svd(bond_matrix);
    S = U * sqrt(Sigma) * V';
    
    % T_{ijk} = sum_m S_{mi} S_{mj} S_{mk}
    T = tensor(zeros(2, 2, 2));
    
    for i=1:2
        for j=1:2
            for k=1:2
                for m=1:2
                    T(i,j,k) = T(i,j,k) + ...
                        S(m,i) * S(m,j) * S(m,k);
                end
            end
        end
    end
    
    
end
