%% Code for AFFGT(2021) Replication

% Linghui Wu

%% Housekeeping
clc;
clear;
close all;


%% Load data
load("../data/AFFGT_data.mat")

params.theta=4;
params.sigma=4;
params.f_d=1;
params.f_u=1;
params.alpha_d=1-data.Input_share_downs(1);
params.alpha_u=1;

params.L_us = 10 * data.L(1) / sum(data.L);
params.L_row = 10 * data.L(2) / sum(data.L);
params.A_d_us = 1;
params.A_u_us = 1;
params.A_d_row = 0.2752;
params.A_u_row = 0.1121;
params.tau_d = 3.0066;
params.tau_u = 2.5971;
params.w_us = 1;


% %% Solve for zero-tariff equilibrium allocation
% 
% % 7 unknown params
% % w_j: x(1)^2
% % M_u_i: x(2)^2, M_u_j: x(3)^2, M_d_i: x(4)^2, M_d_j: x(5)^2
% % T_i: x(6)^2, T_j: x(7)^2
%
% Set initial guesses and optimize
% rng("default");
x0 = unifrnd(0, 1, [7, 1]);  % Vector of unknowns

% Set optimization options
opt1 = optimset("Algorithm", "levenberg-marquardt", "Display", "off", ...
    "MaxFunEvals", 1e4, "MaxIter", 1e4, "TolX", 1e-6, ...
    "TolFun", 1e-6, "Diagnostics", "off");
[x, fval, exitflag1, output] = fsolve(@(x)solve_eqlm(params, x), x0, opt1);
x_squared = x .^ 2;
disp("**********");
disp(exitflag1); % 1: Equation solved. First-order optimality is small.
% Print equilibrium
disp(x_squared);
% Calculates statistics around the zero-tariff equilibrium
% disp(cal_stats(x_squared, params));

%% Solve for optimal tax instruments

% Optimal tariffs
% 9 unknown params
% w_j: x(1)^2
% M_u_i: x(2)^2, M_u_j: x(3)^2, M_d_i: x(4)^2, M_d_j: x(5)^2
% T_i: x(6)^2, T_j: x(7)^2
% t_u_ji: x(8)^2, t_d_ji: x(9)^2;
x0=unifrnd(0,1,[9,1]);
opt2=optimoptions(@fmincon,"Display","off","Algorithm","SQP",...
    "MaxIterations", 2e3,"MaxFunctionEvaluations",2e3,...
    "OptimalityTolerance",1e-6);
% Non-linear constraints
nonlcon_fun=@(x) opt_tariff_cons(x,params);
lb=1e-3*ones(1,9); 
ub=ones(1,9)-1e-3; 
f=@(x) solve_opt_tariff(x,params);
problem=createOptimProblem("fmincon","objective",f,"x0",x0,"lb",lb,...
    "ub",ub,"nonlcon",nonlcon_fun,"options",opt2);
gs=GlobalSearch;
xgs=run(gs,problem);
disp("**********");
display(f(xgs));
display(xgs.^2);
% display(cal_stats(xgs.^2, params));

% Optimal tariffs and domestic subsidy
% 10 unknown params
% w_j: x(1)^2
% M_u_i: x(2)^2, M_u_j: x(3)^2, M_d_i: x(4)^2, M_d_j: x(5)^2
% T_i: x(6)^2, T_j: x(7)^2
% t_u_ji: x(8)^2, t_d_ji: x(9)^2;
% t_u_ii: -x(10)^2
x0=unifrnd(0,1,[10,1]);
opt3=optimoptions(@fmincon,"Display","iter","Algorithm","interior-point",...
    "MaxIterations",1e4,"MaxFunctionEvaluations",1e4,...
    "OptimalityTolerance",1e-6,"ConstraintTolerance",1e-6,...
    "StepTolerance",1e-6);
% Non-linear constraints
nonlcon_fun2=@(x) opt_tax_tar_cons(x,params);
lb=1e-4*ones(1,10);
ub=ones(1,10)-1e-4; 
f2=@(x) solve_opt_tar_sub(x,params);
problem2=createOptimProblem("fmincon","objective",f2,"x0",x0,"lb",lb,...
    "ub",ub,"nonlcon",nonlcon_fun2,"options",opt3);
gs2=GlobalSearch;
xgs2=run(gs2,problem2);
display(f2(xgs2));
display(xgs2.^2);

% Optimal tax policy
% 11 unknown params
% w_j: x(1)^2
% M_u_i: x(2)^2, M_u_j: x(3)^2, M_d_i: x(4)^2, M_d_j: x(5)^2
% T_i: x(6)^2, T_j: x(7)^2
% t_u_ji: x(8)^2, t_d_ji: x(9)^2;
% t_u_ii: -x(10)^2, v_d_ij: x(11)^2;
x0=unifrnd(0,1,[11,1]);
opt4=optimoptions(@fmincon,"Display","final","Algorithm","interior-point",...
    "MaxIterations",1e4,"MaxFunctionEvaluations",1e4,...
    "OptimalityTolerance",1e-6,"ConstraintTolerance",1e-6,...
    "StepTolerance",1e-6);
% Non-linear constraints
nonlcon3=@(x) opt_tax_plc_cons(x,params);
lb=1e-2*ones(1,11);
ub=ones(1,11)-1e-3;
f3=@(x) solve_opt_tax_plc(x,params);
problem3=createOptimProblem("fmincon","objective",f3,"x0",x0,"lb",lb,...
    "ub",ub,"nonlcon",nonlcon3,"options",opt4);
gs3=GlobalSearch;
xgs3=run(gs3,problem3);
display(f3(xgs3));
display(xgs3.^2);

%% Use nested fixed point to solve the optimal tariffs problem

% Solve the zero-tariff equilibrium
% t = zeros(1,8);
% v = zeros(1,8);
% solve_opt_tax(t, v, params);

% Only tariffs: t_u_ji, t_d_ji
% t0 = unifrnd(0, 1, [2,1]);
% f = @(t) solve_opt_tax(t, params);
% option_t1=optimoptions(@fmincon,"Display","off","Algorithm","SQP",...
%     "MaxIterations", 1e4,"OptimalityTolerance",1e-8);
% A = [];
% b = [];
% Aeq = [];
% beq = [];
% lb = 1e-3 * ones(1,2);
% ub = ones(1,2) - 1e-3;
% nonlcon = []; % No further non-linear constraints included
% f=@(t) solve_opt_tax(t,params);
% problem_t1=createOptimProblem("fmincon","objective",f,"x0",t0,"lb",lb,...
%     "ub",ub,"nonlcon",nonlcon,"options",option_t1);
% gs=GlobalSearch;
% xgs=run(gs,problem_t1);
% opt_tax = optimoptions(@fmincon,"Display","off","Algorithm","SQP",...
%         "MaxIterations",1e4,"OptimalityTolerance",1e-8);
% [t,utility,exitflag_tax,output] = fmincon(f,t0,A,b,Aeq,beq,lb,ub,nonlcon,opt_tax);
