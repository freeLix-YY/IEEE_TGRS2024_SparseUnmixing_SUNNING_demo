
clear;
close all;
clc;

set(0,'DefaultFigureWindowStyle','docked')

addpath("DC2");
addpath("Utilities");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate DC2 Abundance Map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% X: 2D abundance maps
% nl: height
% nc: width
% p: number of endmembers
[X,nl,nc,p] = GenerateAbundanceDC2(); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build a synthetic library
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A: synthetic library (bands x number of materials)
% M: GT library (bands x endmember materials)
% supp: index of endmember materials in A
[A, M, supp] = BuildLibraryDC2();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load noise case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load("DC2_10%_Uniform_Noise.mat");
load("DC2_Gauss_SNR=20to35dB.mat");

%  observed spectral vector
Y = M * X + gauss_noise + rand_matrix; % additive model

% create true X wrt the library A
L = size(A,1); % L: total number of bands
n = size(A,2); % n: total number of materials
N = nl*nc; % N: total number of pixels

XT = zeros(n,N); % GT abundance maps, its size followed by A
XT(supp,:) = X; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run sparse unmixing algorithms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath("LRUnSAL_TV");

[LRUnSAL_TV.X_hat,~,~] = LRUnSAL_TV(A, Y, 'MU', 0.05, 'POSITIVITY', 'yes', 'ADDONE', 'no', ...
                                    'LAMBDA_1', 0.0001, 'LAMBDA_TV', 0.01, 'LAMBDA_2',0.01, 'TV_TYPE', 'niso',...
                                    'IM_SIZE', [nl,nc], 'AL_ITERS', 300, 'TRUE_X', XT, 'TOL', 1e-6 , 'VERBOSE', 'no');  

[LRUnSAL_TV.avg_AbunPSNR, LRUnSAL_TV.avg_AbunSSIM, LRUnSAL_TV.avg_AbunSRE] = ...
           Evaluation(LRUnSAL_TV.X_hat, XT, nl, nc, n, supp);

DisplayAbundance(LRUnSAL_TV.X_hat, XT, nl, nc, n, p, supp, 'LRUnSAL-TV');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run SUNNING algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preparation
input.maxIter = 1000; % number of iterations
input.alpha = 20; % step-size
input.hyper_a = 100; % log-cosh
input.sparsity = p; % l0-norm
input.A = A; % spectral library
input.Y = Y; % observed HSI

init.X0 = LRUnSAL_TV.X_hat; % initial point

info.display = true; % display evaluation
info.displayIter = 100; % 100 iterations display once

info.XT = XT; % GT abundance maps
info.nl = nl; % nl: height
info.nc = nc; % nc: width
info.n = n; % n: total number of materials
info.p = p; % p: number of endmembers
info.supp = supp; % supp: index of endmember materials in A

info.normalization = false; % don't need normalization 

[SUNNING.X_hat] = SUNNING(input, init, info);

[SUNNING.avg_AbunPSNR, SUNNING.avg_AbunSSIM, SUNNING.avg_AbunSRE] = ...
           Evaluation(SUNNING.X_hat, XT, nl, nc, n, supp);

DisplayAbundance(SUNNING.X_hat, XT, nl, nc, n, p, supp, 'SUNNING');
        