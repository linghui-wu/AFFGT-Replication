%% Code for AFFGT(2021) Replication

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
% w_i: x(1) ^ 2, w_j: x(2) ^ 2
% M_d_i: x(3) ^ 2, M_d_j: x(4) ^ 2, M_u_i: x(5) ^ 2, M_u_j: x(6) ^ 2
% t_i: x(7) ^ 2, T_j: x(8) ^ 2

% Set initial guesses and optimize
% rng("default");
x0 = unifrnd(0, 1, [8, 1]);  % Vector of unknowns

% Set optimization options
alg0 = "trust-region-dogleg";
alg1 = "trust-region";
alg2 = "levenberg-marquardt";
opt = optimset("Algorithm", alg0, "Display", "iter", "MaxFunEvals", 1e5, ...
    "MaxIter", 1e5, "TolX", 1e-6, "TolFun", 1e-6, "Diagnostics", "on");

x = fsolve(@(x)solve_eq(params, x), x0, opt);
disp(x .^ 2);


%% Check 
% clear;
% clc;
% 
% % Unknown params
% syms w_i w_j M_d_i M_d_j M_u_i M_u_j
% % Known params
% syms theta sigma % Substitution elasticities
% syms f_d f_u % Fixed costs
% syms alpha_d alpha_u % Labor intensity
% syms L_i L_j % Scaled population
% syms A_d_i A_d_j A_u_i A_u_j % Technology
% syms tau_d tau_u % Iceberg trade costs
% syms t_d_ii t_d_ij t_u_ii t_u_ij t_d_jj t_d_ji t_u_jj t_u_ji % Taxes
% syms v_d_ii v_d_ij v_u_ii v_u_ij v_d_jj v_d_ji v_u_jj v_u_ji % Subsidies
% syms T_i T_j % Tax revenues
% 
% % Mark-ups
% mu_u = theta / (theta - 1);
% mu_d = sigma / (sigma - 1);
% 
% % Labor intensity
% alpha_d_bar = 1 / (alpha_d ^ alpha_d * (1 - alpha_d) ^ (1 - alpha_d));
% alpha_u_bar = 1 / (alpha_u ^ alpha_u * (1 - alpha_d) ^ (1 - alpha_u)); 
% 
% 
% %% Equilibrium equations
% 
% % Marginal cost of upstream sector in country i and j
% mc_u_i = alpha_u_bar / A_u_i * w_i ^ alpha_u; % * (P_u_i) ^ 0 
% mc_u_j = alpha_u_bar / A_u_j * w_j ^ alpha_u; % * (P_u_j) ^ 0
% 
% % Price in upstream sectors
% p_u_ij = mu_u * tau_u * mc_u_i / (1 + v_u_ij);
% p_u_ji = mu_u * tau_u * mc_u_j / (1 + v_u_ji);
% p_u_ii = mu_u * mc_u_i / (1 + v_u_ii);
% p_u_jj = mu_u * mc_u_j / (1 + v_u_jj);
% 
% % Price index in upstream sectors
% P_u_ji = M_u_j ^ (1 / (1 - theta)) * (1 + t_u_ji) * p_u_ji;
% P_u_ij = M_u_i ^ (1 / (1 - theta)) * (1 + t_u_ij) * p_u_ij;
% P_u_ii = M_u_i ^ (1 / (1 - theta)) * (1 + t_u_ii) * p_u_ii;
% P_u_jj = M_u_j ^ (1 / (1 - theta)) * (1 + t_u_jj) * p_u_jj;
% 
% % Price index in upstream sector for country i and country j
% P_u_i = (P_u_ii ^ (1 - theta) + P_u_ji ^ (1 - theta)) ^ (1 / (1 - theta));
% P_u_j = (P_u_jj ^ (1 - theta) + P_u_ij ^ (1 - theta)) ^ (1 / (1 - theta));
% 
% % Marginal cost of downstream sectors in country i and j
% mc_d_i = alpha_d_bar / A_d_i * w_i ^ alpha_d * P_u_i ^ (1 - alpha_d);
% mc_d_j = alpha_d_bar / A_d_j * w_j ^ alpha_d * P_u_j ^ (1 - alpha_d);
% 
% % Price in downstream sectors
% p_d_ij = mu_d * tau_d * mc_d_i / (1 + v_d_ij);
% p_d_ji = mu_d * tau_d * mc_d_j / (1 + v_d_ji);
% p_d_ii = mu_d * mc_d_i / (1 + v_d_ii);
% p_d_jj = mu_d * mc_d_i / (1 + v_d_jj);
% 
% % Price index in downstream sectors
% P_d_ji = M_u_j ^ (1 / (1 - sigma)) * (1 + t_d_ji) * p_d_ji;
% P_d_ij = M_u_i ^ (1 / (1 - sigma)) * (1 + t_d_ij) * p_d_ij;
% P_d_ii = M_u_i ^ (1 / (1 - sigma)) * (1 + t_d_ii) * p_d_ii;
% P_d_jj = M_u_i ^ (1 / (1 - sigma)) * (1 + t_d_jj) * p_d_jj;
% 
% % Price index in downstream sectors in country i and j
% P_d_i = (P_d_ii ^ (1 - sigma) + P_d_ji ^ (1 - sigma)) ^ (1 / (1 - sigma));
% P_d_j = (P_d_jj ^ (1 - sigma) + P_d_ij ^ (1 - sigma)) ^ (1 / (1 - sigma));
% 
% % Free entry
% y_d_i = (sigma - 1) * f_d;
% y_d_j = (sigma - 1) * f_d;
% y_u_i = (theta - 1) * f_u;
% y_u_j = (theta - 1) * f_u;
% 
% % Labor demand for upstream and downstream sectors
% l_d_i = alpha_d * mc_d_i * (f_d + y_d_i) / w_i;
% l_d_j = alpha_d * mc_d_j * (f_d + y_d_j) / w_j;
% l_u_i = (f_u + y_u_i) / A_u_i;
% l_u_j = (f_u + y_u_j) / A_u_j;
% 
% % Quantities in upstream sectors
% Q_u_ji = (1 - alpha_d) * mc_d_i * (f_d + y_d_i) / P_u_i * (P_u_ji / P_u_i) ^ (-theta);
% Q_u_jj = (1 - alpha_d) * mc_d_j * (f_d + y_d_j) / P_u_j * (P_u_jj / P_u_j) ^ (-theta);
% Q_u_ij = (1 - alpha_d) * mc_d_j * (f_d + y_d_j) / P_u_j * (P_u_ij / P_u_j) ^ (-theta);
% Q_u_ii = (1 - alpha_d) * mc_d_i * (f_d + y_d_i) / P_u_i * (P_u_ii / P_u_i) ^ (-theta);
% 
% % Output in upstream sectors
% x_ji = Q_u_ji * ((1 + t_u_ji) * p_u_ji / P_u_ji) ^ (-theta);
% x_jj = Q_u_jj * ((1 + t_u_jj) * p_u_jj / P_u_jj) ^ (-theta);
% x_ij = Q_u_ij * ((1 + t_u_ij) * p_u_ij / P_u_ij) ^ (-theta);
% x_ii = Q_u_ii * ((1 + t_u_ii) * p_u_ii / P_u_ii) ^ (-theta);
% 
% % Output in downstream sectors
% c_ji = (w_i * L_i + T_i) / P_d_i ^ (1 - sigma) * ((1 + t_d_ji) * p_d_ji) ^ (-sigma);
% c_jj = (w_j * L_j + T_j) / P_d_j ^ (1 - sigma) * ((1 + t_d_jj) * p_d_jj) ^ (-sigma);
% c_ij = (w_j * L_j + T_j) / P_d_j ^ (1 - sigma) * ((1 + t_d_ij) * p_d_ij) ^ (-sigma);
% c_ii = (w_i * L_i + T_i) / P_d_i ^ (1 - sigma) * ((1 + t_d_ii) * p_d_ii) ^ (-sigma);
% 
% 
% %% Equilibrium conditions
% 
% % Labor market clearing
% LMC_i = L_i - M_u_i * l_d_i - M_d_i * l_u_i;
% LMC_j = L_j - M_u_j * l_d_j - M_d_j * l_u_j;
% 
% % Goods market clearing
% GMC_u_i = y_u_i - M_u_i * x_ii - M_u_j * tau_u * x_ij;
% GMC_u_j = y_u_j - M_u_j * x_jj - M_u_i * tau_u * x_ji;
% GMC_d_i = y_d_i - c_ii - tau_d * c_ji;
% GMC_d_j = y_d_j - c_jj - tau_d * c_ij;
% 
% eq = [LMC_i, LMC_j, ...
%     GMC_u_i, GMC_u_j, GMC_d_i, GMC_d_j, ...
%     T_i, T_j];

