function eq = solve_eq(params, x)

% disp(params);

%% Set params
% Notations: d for downstream, u for upstream; i for us, j for row.

% Substitution elasticities
theta = params.theta;
sigma = params.sigma;

% Mark-ups
mu_u = theta / (theta - 1);
mu_d = sigma / (sigma - 1);

% Fixed costs (country independent)
f_d = params.f_d;
f_u = params.f_u;

% Labor intensity
alpha_d = params.alpha_d;
alpha_u = params.alpha_u;

alpha_d_bar = 1 / (alpha_d ^ alpha_d * (1 - alpha_d) ^ (1 - alpha_d));
alpha_u_bar = 1 / (alpha_u ^ alpha_u * (1 - alpha_u) ^ (1 - alpha_u)); 

% Population
L_i = params.L_us;
L_j = params.L_row;

% Technology
A_d_i = params.A_d_us;
A_d_j = params.A_d_row;
A_u_i = params.A_u_us;
A_u_j = params.A_u_row;

% Iceberg trade costs
tau_d = params.tau_d;
tau_u = params.tau_u;

% Wages: normalize US wage to 1
% w_i = params.w_us;

% Taxes (t_s_ij rep'ts tax by country j on goods from country i sector s)
% Assumes no tariffs on US goods in baseline
t_d_ii = 0;  % Local tax on downstream goods in country i
t_d_ij = 0;
t_d_jj = 0;  % Local tax on downstream goods in country j
t_d_ji = 0;
t_u_ii = 0;  % Local tax on upstream goods in country i
t_u_ij = 0;
t_u_jj = 0;  % Local tax on upstream goods in country 
t_u_ji = 0;

% Subsidies (v_s_ji denotes subsidies to goods from country j sector s to country i)
v_d_ii = 0;
v_d_ij = 0;
v_d_jj = 0;
v_d_ji = 0;
v_u_ii = 0;
v_u_ij = 0;
v_u_jj = 0;
v_u_ji = 0;


%% Equilibrium equations

% Marginal cost of upstream sector in country i and j
mc_u_i = alpha_u_bar / A_u_i * w_i ^ alpha_u; % * (P_u_i) ^ 0 
mc_u_j = alpha_u_bar / A_u_j * w_j ^ alpha_u; % * (P_u_j) ^ 0

% Price in upstream sectors
p_u_ij = mu_u * tau_u * mc_u_i / (1 + v_u_ij);
p_u_ji = mu_u * tau_u * mc_u_j / (1 + v_u_ji);
p_u_ii = mu_u * mc_u_i / (1 + v_u_ii);
p_u_jj = mu_u * mc_u_j / (1 + v_u_jj);

% Price index in upstream sectors
P_u_ji = M_u_j ^ (1 / (1 - theta)) * (1 + t_u_ji) * p_u_ji;
P_u_ij = M_u_i ^ (1 / (1 - theta)) * (1 + t_u_ij) * p_u_ij;
P_u_ii = M_u_i ^ (1 / (1 - theta)) * (1 + t_u_ii) * p_u_ii;
P_u_jj = M_u_j ^ (1 / (1 - theta)) * (1 + t_u_jj) * p_u_jj;

% Price index in upstream sector for country i and country j
P_u_i = (P_u_ii ^ (1 - theta) + P_u_ji ^ (1 - theta)) ^ (1 / (1 - theta));
P_u_j = (P_u_jj ^ (1 - theta) + P_u_ij ^ (1 - theta)) ^ (1 / (1 - theta));

% Marginal cost of downstream sectors in country i and j
mc_d_i = alpha_d_bar / A_d_i * w_i ^ alpha_d * P_u_i ^ (1 - alpha_d);
mc_d_j = alpha_d_bar / A_d_j * w_j ^ alpha_d * P_u_j ^ (1 - alpha_d);

% Price in downstream sectors
p_d_ij = mu_d * tau_d * mc_d_i / (1 + v_d_ij);
p_d_ji = mu_d * tau_d * mc_d_j / (1 + v_d_ji);
p_d_ii = mu_d * mc_d_i / (1 + v_d_ii);
p_d_jj = mu_d * mc_d_j / (1 + v_d_jj);

% Price index in downstream sectors
P_d_ji = M_d_j ^ (1 / (1 - sigma)) * (1 + t_d_ji) * p_d_ji;
P_d_ij = M_d_i ^ (1 / (1 - sigma)) * (1 + t_d_ij) * p_d_ij;
P_d_ii = M_d_i ^ (1 / (1 - sigma)) * (1 + t_d_ii) * p_d_ii;
P_d_jj = M_d_j ^ (1 / (1 - sigma)) * (1 + t_d_jj) * p_d_jj;

% Price index in downstream sectors in country i and j
P_d_i = (P_d_ii ^ (1 - sigma) + P_d_ji ^ (1 - sigma)) ^ (1 / (1 - sigma));
P_d_j = (P_d_jj ^ (1 - sigma) + P_d_ij ^ (1 - sigma)) ^ (1 / (1 - sigma));

% Free entry
y_d_i = (sigma - 1) * f_d; y_d_j = y_d_i;
y_u_i = (theta - 1) * f_u; y_u_j = y_u_i;

% Labor demand for upstream and downstream sectors
l_d_i = alpha_d * mc_d_i * (f_d + y_d_i) / w_i;
l_d_j = alpha_d * mc_d_j * (f_d + y_d_j) / w_j;
l_u_i = (f_u + y_u_i) / A_u_i;
l_u_j = (f_u + y_u_j) / A_u_j;

% Quantities in upstream sectors
Q_u_ji = (1 - alpha_d) * mc_d_i * (f_d + y_d_i) / P_u_i * (P_u_ji / P_u_i) ^ (-theta);
Q_u_jj = (1 - alpha_d) * mc_d_j * (f_d + y_d_j) / P_u_j * (P_u_jj / P_u_j) ^ (-theta);
Q_u_ij = (1 - alpha_d) * mc_d_j * (f_d + y_d_j) / P_u_j * (P_u_ij / P_u_j) ^ (-theta);
Q_u_ii = (1 - alpha_d) * mc_d_i * (f_d + y_d_i) / P_u_i * (P_u_ii / P_u_i) ^ (-theta);

% Output in upstream sectors
x_ji = Q_u_ji * ((1 + t_u_ji) * p_u_ji / P_u_ji) ^ (-theta);
x_jj = Q_u_jj * ((1 + t_u_jj) * p_u_jj / P_u_jj) ^ (-theta);
x_ij = Q_u_ij * ((1 + t_u_ij) * p_u_ij / P_u_ij) ^ (-theta);
x_ii = Q_u_ii * ((1 + t_u_ii) * p_u_ii / P_u_ii) ^ (-theta);

% Output in downstream sectors 
% No tax revenues
c_ji = (w_i * L_i + T_i) * ((1 + t_d_ji) * p_d_ji) ^ (-sigma) / P_d_i ^ (1 - sigma);
c_jj = (w_j * L_j + T_j) * ((1 + t_d_jj) * p_d_jj) ^ (-sigma) / P_d_j ^ (1 - sigma);
c_ij = (w_j * L_j + T_j) * ((1 + t_d_ij) * p_d_ij) ^ (-sigma) / P_d_j ^ (1 - sigma);
c_ii = (w_i * L_i + T_i) * ((1 + t_d_ii) * p_d_ii) ^ (-sigma) / P_d_i ^ (1 - sigma);


%% Equilibrium constraints

% Labor market clearing
LMC_i = L_i - M_d_i * l_d_i - M_u_i * l_u_i;
LMC_j = L_j - M_d_j * l_d_j - M_u_j * l_u_j;

% Goods market clearinge
GMC_d_i = y_d_i - c_ii - tau_d * c_ij;
GMC_d_j = y_d_j - c_jj - tau_d * c_ji;
GMC_u_i = y_u_i - M_d_i * x_ii - M_d_j * tau_u * x_ij;
GMC_u_j = y_u_j - M_d_j * x_jj - M_d_i * tau_u * x_ji;

% Buget balance in country i and j
% To be completed if not zero-tariff equilibrium
BB_i = T_i;
BB_j = T_j;

%% Function output
eq = [LMC_i, LMC_j, ...
    GMC_u_i, GMC_u_j, GMC_d_i, GMC_d_j, ...
    BB_i, BB_j];

% disp(eq);

end