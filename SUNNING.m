%% The MATLAB implementation of SUNNING
% 
% If you use this code, please cite the following paper:
% 
%    @article{chan2024sparse,
%      title={Sparse Unmixing in the Presence of Mixed Noise Using $\ell$_0-Norm Constraint and Log-Cosh Loss},
%      author={Chan, Yiu Yu and Li, Xiao-Peng and Mai, Jiajie and Leung, Chi-Sing and So, Hing Cheung},
%      journal={IEEE Transactions on Geoscience and Remote Sensing},
%      year={2024},
%      publisher={IEEE}
%    }

% This code solves the L0-norm sparse unmixing optimization problem:
% 
%  min  sum(log(cosh(a * (A * X - Y))), "all") / a
%   X
% subject to X(:,j) >= 0, 1^T*X(:,j) = 1, ||X(:,j)||_0 <= S, for j = [1, N-pixels]

%% Input:
% A: spectral library (L bands x n materials)
% Y: observed spectral vector (L bands x N-pixels)
% S: l0-norm constraint, limited to the number of existing endmembers
% a: hyper-parameters to log-cosh loss 
% alpha: step size to gradient descent
% maxIter: maximum number of iterations

%% Ouput:
% X: abundance matrix (n materials x N-pixels)

%% (If needed) Evaluation:
% XT: ground truth abundance matrix (p endmembers x N-pixels)
% nl: height
% nc: width 
% n: total number of materials
% p: number of endmembers
% supp: index of endmember materials in A
% normalization: used for Samson, JasperRidge, Urban datasets
% display: true = display evaluation
% displayIter: i-th iterations display once

function [X] = SUNNING(input, init, info)

% assign variables
[A, Y, S, a, alpha, maxIter] = deal(input.A, input.Y, input.sparsity, input.hyper_a, input.alpha, input.maxIter);
clearvars input

if isempty(init.X0)
    % We use any of sparse unmixing algorithms as initial X^0.
    % For example, we use RHUIDR_HTV or LRUnSAL_TV as coarse solution.
    % Then, SUNNING can be implemented as fine stage method.

    % The MATLAB implementation of RHUIDR_HTV:
    % links: https://github.com/MDI-TokyoTech/Towards_Robust_Hyperspectral_
    % Unmixing_Mixed_Noise_Modeling_and_Image-Domain_Regularization

    % (Hyper)-parameters:
    % MU : Lagrangian weight
    % POSITIVITY : Non-negativity constraint to X
    % ADDONE : Sum-to-one constraint to X
    % LAMBDA_1 : Regularization to nuclear norm of X
    % LAMBDA_TV : Regularization to total variation of X
    % LAMBDA_2 : Regularization to L1-sparse noise of S

    [X,~,~] = LRUnSAL_TV(A, Y, 'MU', 0.05, 'POSITIVITY', 'yes', 'ADDONE', 'no', ...
                        'LAMBDA_1', 0.0001, 'LAMBDA_TV', 0.001, 'LAMBDA_2', 0.01, 'TV_TYPE', 'niso',...
                        'IM_SIZE', [info.nl,info.nc], 'AL_ITERS', 1000, 'TOL', 1e-6 , 'VERBOSE', 'no');  

else
    X = init.X0; % initial point
    clearvars init
end

% L-smooth bound
L = a * max(eig(A'*A));
step_size = alpha/L;

% function handle
G_gradient = @(X) A' * tanh(a * (A * X - Y));
gradient_descent = @(X) X - step_size * G_gradient(X);

% temporary variables
zero_x = zeros(info.n,1);

for t = 1:maxIter

    % two-step alternating scheme: gradient descent + projection
    X = projection(gradient_descent(X), S, zero_x); 

    % (if needed) performance evaluation 
    if (info.display == true) 
        if (mod(t, info.displayIter) == 0) || t == 1
            if info.normalization == false
                [avg_AbunPSNR, avg_AbunSSIM, avg_AbunSRE] = Evaluation(X, info.XT, info.nl, info.nc, info.n, info.supp);
            else
                X_avg = [];
                for i = 1:info.p
                    X_avg = [X_avg; 
                             mean(X(info.A_sub_index(i,1):info.A_sub_index(i,2),:),1)];
                end
    
                % normalization to fit sum-to-one
                X_avg = X_avg ./ repmat(sum(X_avg, 1), [info.p 1]);
    
                [avg_AbunPSNR, avg_AbunSSIM, avg_AbunSRE] = Evaluation(X_avg, info.XT, info.nl, info.nc, info.p, info.supp);
            end
            fprintf("\nAt %d-th iteration :\n", t);
            fprintf("SRE = %2.3f, PSNR = %2.3f, SSIM = %2.3f\n", avg_AbunSRE, avg_AbunPSNR, avg_AbunSSIM);
        end
    end
end
end


function [Z] = projection(Z, S, pre_zero_x_vec)
    
    [Z_sort, Mat_sort_idx] = sort(Z, 1, 'descend'); % sort each column vector in descending order

    % only care the sorted elements before S-th position or at S-th position
    Z_sub_sort = Z_sort(1:S,:); 

    % check z{s+1} <= mu{s}/2
    parfor i = 1:size(Z,2)
        Z(:,i) = searching(Z_sub_sort(:,i), S, Mat_sort_idx(:,i), pre_zero_x_vec);
    end
end


function x = searching(z_sort, S, recover_idx, x)
    % initial mu{1}
    mu = 2*(z_sort(1)-1);
  
    % check k = 1, z{k+1} <= mu{k}/2
    if (z_sort(2) <= mu/2)
        x(1) = z_sort(1) - mu/2;
        x(recover_idx,1) = x; 
        return
    end

    % check k = [2,...,S-1] 
    for k = 2: S-1
        mu = ((k-1) * mu + 2 * z_sort(k)) / k; 
        if (z_sort(k + 1) <= mu/2) 
            % Case 1: largest possible sparsity k
            % sorted_x {1,...,k} keeps k-th largest elements
            % input is pre-ready for x{k+1,...,n} sets to zero, more efficient
            x(1:k) = z_sort(1:k) - mu/2;
            x(recover_idx,1) = x; % sorted solution x push back to orginal order
            return
        end
    end
    
    % Case 2: largest sparsity S, i.e. consider z{s} > mu{s-1}/2
    x(1:S) = z_sort(1:S) - (((S-1) * mu + 2 * z_sort(S)) / S) / 2;
    x(recover_idx,1) = x; % sorted solution x push back to orginal order
   
end