%% Code for AFFGT(2021) Replication

% Linghui Wu

%% Housekeeping
clear;
close all;


%% Load data
load("../AFFGT_data.mat")

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
% % Seven unknown params
% % w_j: x(1)^2
% % M_u_i: x(2)^2, M_u_j: x(3)^2, M_d_i: x(4)^2, M_d_j: x(5)^2
% % T_i: x(6)^2, T_j: x(7)^2
% 
% % Set initial guesses and optimize
% % rng("default");
% x0 = unifrnd(0, 1, [7, 1]);  % Vector of unknowns
% 
% % Set optimization options
% opt1 = optimset("Algorithm", "levenberg-marquardt", "Display", "off", ...
%     "MaxFunEvals", 1e4, "MaxIter", 1e4, "TolX", 1e-6, ...
%     "TolFun", 1e-6, "Diagnostics", "off");
% [x, fval, exitflag1, output] = fsolve(@(x)solve_eqlm(params, x), x0, opt1);
% x_squared = x .^ 2;
% disp("**********");
% disp(exitflag1); % 1: Equation solved. First-order optimality is small.
% 
% % Print equilibrium
% disp(["Wages", 1, x_squared(1)]);
% disp(["Measure of firms", reshape(x_squared(2:5), [1, 4])]);
% disp(["Tax revenues", reshape(x_squared(6:7), [1, 2])])
% 
% % Calculates statistics around the zero-tariff equilibrium
% disp(cal_stats(x_squared, params));


%% Solve for optimal tax instruments

% Optimal tariffs
% Nine unknown params
% w_j: x(1)^2
% M_u_i: x(2)^2, M_u_j: x(3)^2, M_d_i: x(4)^2, M_d_j: x(5)^2
% T_i: x(6)^2, T_j: x(7)^2
% t_u_ji: x(8)^2, t_d_ji: x(9)^2;

fval=0;
% exitflag2="-1";
cons=zeros(1,7);
while (fval>-0.0317) || (exitflag2<0) || any(cons<=1e-3)
    x0=unifrnd(0,1,[9,1]);
    opt2=optimoptions(@fmincon,"Display","off","MaxIter", 1e4);
    % Non-linear constraints
    nonlcon=@(x) opt_tariff_cons(x,params);
    lb=1e-3*ones(1,9);
    ub=ones(1,9)-1e-3;
    f=@(x) solve_opt_tariff(x,params);
    % No other constraints: A,b,Aeq,beq=[]
    [x,fval,exitflag2]=fmincon(f,x0,[],[],[],[],lb,ub,nonlcon,opt2);
    cons=[x(1),x(2),x(3),x(4),x(5),x(8),x(9)].^2;
    disp(cons);
    disp(fval);
end
disp(x.^2);
