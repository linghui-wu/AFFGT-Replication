%% Code for AFFGT(2021) Replication

% Linghui Wu

%% Housekeeping
clear;
close all;


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
% w_i: x(1) ^ 2, w_j: x(2) ^ 2
% M_u_i: x(3) ^ 2, M_u_j: x(4) ^ 2, M_d_i: x(5) ^ 2, M_d_j: x(6) ^ 2
% T_i: x(7) ^ 2, T_j: x(8) ^ 2
% c_ij: x(7) ^ 2, c_ii: x(8) ^ 2, c_ji: x(9) ^ 2, c_jj: x(10) ^ 2

% w_i: x(1) ^ 2
% M_u_i: x(2) ^ 2, M_u_j: x(3) ^ 2, M_d_i: x(4) ^ 2, M_d_j: x(5) ^ 2
% T_i: x(6) ^ 2, T_j: x(7) ^ 2

% Set initial guesses and optimize
% rng("default");
x0 = unifrnd(0, 1, [8, 1]);  % Vector of unknowns

% Set optimization options
alg0 = "trust-region-dogleg";
alg1 = "trust-region";
alg2 = "levenberg-marquardt";
opt = optimset("Algorithm", alg0, "Display", "off", "MaxFunEvals", 3e5, ...
    "MaxIter", 3e5, "TolX", 1e-6, "TolFun", 1e-6, "Diagnostics", "on");

x = fsolve(@(x)solve_eq(params, x), x0, opt);
x_squared = x .^ 2;
disp(x_squared);

% Eight params: M's converge
% w_i: x(1) ^ 2, w_j: x(2) ^ 2
% M_u_i: x(3) ^ 2, M_u_j: x(4) ^ 2, M_d_i: x(5) ^ 2, M_d_j: x(6) ^ 2
% t_i: x(7) ^ 2, T_j: x(8) ^ 2

