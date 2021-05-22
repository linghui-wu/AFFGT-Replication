%% Code for AFFGT(2021) Replication

% Linghui Wu

%% Housekeeping
clear;
close all;


%% Load data
% load("../AFFGT_data.mat")

params.theta = 4;
params.sigma = 4;
params.f_d = 1;
params.f_u = 1;
params.alpha_d = 0.5483; % 1 - data.Input_share_downs(1);
params.alpha_u = 1;

params.L_us = 0.4531; % 10 * data.L(1) / sum(data.L);
params.L_row = 9.5469; % 10 * data.L(2) / sum(data.L);
params.A_d_us = 1;
params.A_u_us = 1;
params.A_d_row = 0.2752;
params.A_u_row = 0.1121;
params.tau_d = 3.0066;
params.tau_u = 2.5971;
params.w_us = 1;


% %% Solve for zero-tariff equilibrium allocation
% 
% % Seven unknown params
% % w_j: x(1) ^ 2
% % % M_u_i: x(2) ^ 2, M_u_j: x(3) ^ 2, M_d_i: x(4) ^ 2, M_d_j: x(5) ^ 2
% % T_i: x(6) ^ 2, T_j: x(7) ^ 2
% 
% % Set initial guesses and optimize
% % rng("default");
% x0 = unifrnd(0, 1, [7, 1]);  % Vector of unknowns
% 
% % Set optimization options
% alg0 = "trust-region-dogleg";
% alg1 = "trust-region";
% alg2 = "levenberg-marquardt";
% opt = optimset("Algorithm", alg2, "Display", "off", "MaxFunEvals", 1e4, ...
%     "MaxIter", 1e4, "TolX", 1e-6, "TolFun", 1e-6, "Diagnostics", "off");
% 
% [x, fval, exitflag, output] = fsolve(@(x)solve_eqlm(params, x), x0, opt);
% x_squared = x .^ 2;
% 
% disp("**********");
% disp(exitflag); % 1: Equation solved. First-order optimality is small.
% 
% 
% % Print equilibrium
% disp(["Wages", 1, x_squared(1)]);
% disp(["Measure of firms", reshape(x_squared(2:5), [1, 4])]);
% disp(["Tax revenues", reshape(x_squared(6:7), [1, 2])])
% 
% % Calculates statistics around the zero-tariff equilibrium
% disp(cal_stats(x_squared, params));


%% Solve for optimal tax instruments

% % Example
% fun = @(x)3*x(1)^2 + 2*x(1)*x(2) + x(2)^2 - 4*x(1) + 5*x(2);
% x0 = unifrnd(0, 1, [2, 1]);
% options = optimoptions(@fminunc,"Display","iter","Algorithm","quasi-newton");
% [x,fval] = fminunc(fun,x0,options);

x0=unifrnd(0,1,[9,1]);
options=optimoptions(@fminunc,"Display","iter",...
    "Algorithm","quasi-newton");
f=@(x) solve_opt_tariff(params,x)
[x,fval,exit]=fminunc(f,x0,params);


% function J = computeCostMulti(X, y, theta);
%   J = 1/(2*size(y,1)) * sum(((X * theta)-y).^2);
% end;
 
% f = @(x,y,z) function(tSet,preco,theta)
% [optTheta, fVal, exit] = fmiunc(f,tSet,preco,theta,options);



