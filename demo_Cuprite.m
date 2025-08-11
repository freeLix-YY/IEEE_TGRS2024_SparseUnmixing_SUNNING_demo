
clear;
close all;
clc;

set(0,'DefaultFigureWindowStyle','docked')

addpath("Cuprite");
addpath("Utilities");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load Cuprite HSI
% download source : https://github.com/ricardoborsoi/MUA_SparseUnmixing/tree/master/real_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load("cuprite_ref.mat");

Y = x; % HSI

nl = Lines; % nl: height 
nc = Columns; % nc: width

clearvars x Lines Columns wavlen

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build a mineral library
% download source : https://github.com/ricardoborsoi/MUA_SparseUnmixing/tree/master/real_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A: mineral library (bands x number of materials)
A = BuildLibraryCuprite();
A = A(BANDS,:); % remove bands of distortion

clearvars BANDS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load noise case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load("Cuprite_10%_Uniform_Noise.mat");

%  observed spectral vector
Y = Y + rand_matrix; % additive model

% create true X wrt the library A
L = size(A,1); % L: total number of bands
n = size(A,2); % n: total number of materials
N = nl*nc; % N: total number of pixels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run sparse unmixing algorithms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath("LRUnSAL_TV");

[LRUnSAL_TV.X_hat,~,~] = LRUnSAL_TV(A, Y, 'MU', 0.05, 'POSITIVITY', 'yes', 'ADDONE', 'no', ...
                                    'LAMBDA_1', 0.001, 'LAMBDA_TV', 0.001, 'LAMBDA_2',0.01, 'TV_TYPE', 'niso',...
                                    'IM_SIZE', [nl,nc], 'AL_ITERS', 2000, 'TOL', 1e-6 , 'VERBOSE', 'yes');  

DisplayAbundanceCuprite(LRUnSAL_TV.X_hat, nl, nc, n, 'LRUnSAL-TV');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run SUNNING algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preparation
input.maxIter = 15000; % number of iterations
input.alpha = 20; % step-size
input.hyper_a = 100; % log-cosh
input.sparsity = 12; % l0-norm
input.A = A; % spectral library
input.Y = Y; % observed HSI

init.X0 = LRUnSAL_TV.X_hat; % initial point

info.display = false; % display evaluation

info.nl = nl; % nl: height
info.nc = nc; % nc: width
info.n = n; % n: total number of materials

[SUNNING.X_hat] = SUNNING(input, init, info);

DisplayAbundanceCuprite(SUNNING.X_hat, nl, nc, n, 'SUNNING');        