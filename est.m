%% Code for AFFGT 2021 Replication

% Linghui Wu

%% Housekeeping
clear;
close all;
clc;


%% Load data
load("../AFFGT_data.mat")

params.theta = 4;
params.sigma = 4;
params.f_d = 1;
params.f_u = 1;
params.alpha_d = 1 - data.Input_share_downs(1);
params.alpha_u = 1;

params.L_us = 10 * data.L(1) / sum(data.L);
params.L_row = 10 * data.L(2) / sum(data.L);
params.A_d_us = 1;
params.A_u_us = 1;
params.A_d_row = 0.2752;
params.A_u_row = 0.1121;
params.tau_d = 3.0066;
params.tau_u = 2.5971;
params.w_us = 1;


%% Solve for equilibrium allocation

% Unknowns
% w_j: log(x(1))
% M_d_i: log(x(2)), M_d_j: log(x(3)), M_u_i: log(x(4)), M_u_j: log(x(5))

% Set initial guesses and optimize
% rng("default");
x0 = unifrnd(0, 1, [5, 1]);  % Vector of unknowns

% Set optimization options
alg0 = "trust-region-dogleg";
alg1 = "trust-region";
alg2 = "levenberg-marquardt";
opt = optimset("Algorithm", alg0, "Display", "iter", "MaxFunEvals", 2000, ...
    "Maxiter", 2000, "TolX", 1e-6, "TolFun", 1e-6, "Diagnostics", "on");

x = fsolve(@(x)solve_eq(params, x), x0, opt);
disp(exp(x));




