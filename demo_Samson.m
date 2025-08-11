
clear;
close all;
clc;

set(0,'DefaultFigureWindowStyle','docked')

addpath("Samson");
addpath("Utilities");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load Samson HSI
% download source : https://github.com/ricardoborsoi/MUA_SparseUnmixing/tree/master/real_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load("samson_1.mat");

Y = V; % HSI

L = nBand; % L: total number of bands
nl = nRow; % nl: height 
nc = nCol; % nc: width
N = nl*nc; % N: total number of pixels

clearvars V nBand nRow nCol

p = 3; % p: number of endmembers


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load abundance map labels
% download source 1: https://lesun.weebly.com/hyperspectral-data-set.html
% download source 2: https://blog.csdn.net/x5675602/article/details/89185854

% For more information, see the following:
% Hyperspectral Unmixing: Ground Truth Labeling, Datasets, Benchmark Performances and Survey
% links: https://arxiv.org/abs/1708.05125
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load("Samson_GT.mat"); % M: GT library, XT: 2D abundance maps

%% Remarks:
%% You can use XT for a rough-only estimation in SRE, SSIM, PSNR to best tune 
%% the model parameters, but they are not real due to endmember variability
%% We only use them for visual comparison in the paper.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the spectral library
% download source : https://github.com/ricardoborsoi/MUA_SparseUnmixing/tree/master/real_data

% For more information, see the following:
% A Fast Multiscale Spatial Regularization for Sparse Hyperspectral Unmixing
% (Full-version) links: https://arxiv.org/abs/1712.01770
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load("spectral_library_samson.mat"); % spectra library: A

n = size(A,2); % n: total number of materials
nLib1 = size(lib1,2); % number of signatures in endmember 1 (Soil) library
nLib2 = size(lib2,2); % number of signatures in endmember 2 (Tree) library
nLib3 = size(lib3,2); % number of signatures in endmember 3 (Water) library

clearvars lib1 lib2 lib3 material_names

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load noise case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load("Samson_10%_Uniform_Noise.mat");

%  observed spectral vector
Y = Y + rand_matrix; % additive model


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run sparse unmixing algorithms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath("LRUnSAL_TV");

[LRUnSAL_TV.X_hat,~,~] = LRUnSAL_TV(A, Y, 'MU', 0.05, 'POSITIVITY', 'yes', 'ADDONE', 'no', ...
                                   'LAMBDA_1', 0.1, 'LAMBDA_TV', 0.001, 'LAMBDA_2', 0.01, 'TV_TYPE', 'niso',...
                                   'IM_SIZE', [nl,nc], 'AL_ITERS', 2000, 'TOL', 1e-6 ,'VERBOSE', 'yes');

% taking average in each endmember library
LRUnSAL_TV.X_avg = [mean(LRUnSAL_TV.X_hat(1:nLib1,:),1);...
                    mean(LRUnSAL_TV.X_hat((nLib1+1):(nLib1+nLib2),:),1); ...
                    mean(LRUnSAL_TV.X_hat((nLib1+nLib2+1):(nLib1+nLib2+nLib3),:),1)];

% normalization to fit sum-to-one
LRUnSAL_TV.X_avg = LRUnSAL_TV.X_avg ./ repmat(sum(LRUnSAL_TV.X_avg, 1), [p 1]);

[LRUnSAL_TV.avg_AbunPSNR, LRUnSAL_TV.avg_AbunSSIM, LRUnSAL_TV.avg_AbunSRE] = ...
           Evaluation(LRUnSAL_TV.X_avg, XT, nl, nc, p, 1:p);

DisplayAbundance(LRUnSAL_TV.X_avg, XT, nl, nc, p, p, 1:p, 'LRUnSAL-TV');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run SUNNING algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preparation
input.maxIter = 2000; % 2000 for demo, 1e4 < maxIter < 1e5 for experiments
input.alpha = 15; % step-size
input.hyper_a = 100; % log-cosh
input.sparsity = 10; % l0-norm, even S > p, it eventually reduces to S = K
input.A = A; % spectral library
input.Y = Y; % observed HSI

init.X0 = LRUnSAL_TV.X_hat; % initial point

info.display = true; % display evaluation
info.displayIter = 1000; % 1000 iterations display once

info.XT = XT; % GT abundance maps
info.nl = nl; % nl: height
info.nc = nc; % nc: width
info.n = n; % n: total number of materials
info.supp = 1:p; % supp: index of endmember materials in A
info.p = p; % p: number of endmembers

info.normalization = true; % need normalization 
info.A_sub_index = [1, nLib1;
                    nLib1+1, nLib1+nLib2;
                    nLib1+nLib2+1, nLib1+nLib2+nLib3]; % index of sub-class of A

[SUNNING.X_hat] = SUNNING(input, init, info);

% taking average in each endmember library
SUNNING.X_avg = [mean(SUNNING.X_hat(1:nLib1,:),1);...
                 mean(SUNNING.X_hat((nLib1+1):(nLib1+nLib2),:),1); ...
                 mean(SUNNING.X_hat((nLib1+nLib2+1):(nLib1+nLib2+nLib3),:),1)];

% normalization to fit sum-to-one
SUNNING.X_avg = SUNNING.X_avg ./ repmat(sum(SUNNING.X_avg, 1), [p 1]);

[SUNNING.avg_AbunPSNR, SUNNING.avg_AbunSSIM, SUNNING.avg_AbunSRE] = ...
           Evaluation(SUNNING.X_avg, XT, nl, nc, p, 1:p);

DisplayAbundance(SUNNING.X_avg, XT, nl, nc, p, p, 1:p, 'SUNNING');