% Computes and compares sum stats around equilibrium
% Retrieve values of optimized params
x = xgs3 .^ 2;
[w_j] = deal(x(1));
[M_u_i, M_u_j, M_d_i, M_d_j] = deal(x(2), x(3), x(4), x(5));
[T_i, T_j] = deal(x(6), x(7));

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
w_i = params.w_us;
% Taxes (t_s_ij rep'ts tax by country j on goods from country i sector s)
% Assumes no tariffs on US goods in baseline
t_d_ii = 0; t_d_ij = 0; t_d_jj = 0; t_d_ji = x(9);
t_u_ii = -x(10); t_u_ij = 0; t_u_jj = 0; t_u_ji = x(8);
% Subsidies (v_s_ji denotes subsidies to goods from country j sector s to country i)
v_d_ii = 0; v_d_ij = 0; v_d_jj = 0; v_d_ji = 0;
v_u_ii = 0; v_u_ij = 0; v_u_jj = 0; v_u_ji = 0;
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
c_ji = (w_i * L_i + T_i) * ((1 + t_d_ji) * p_d_ji) ^ (-sigma) / P_d_i ^ (1 - sigma);
c_jj = (w_j * L_j + T_j) * ((1 + t_d_jj) * p_d_jj) ^ (-sigma) / P_d_j ^ (1 - sigma);
c_ij = (w_j * L_j + T_j) * ((1 + t_d_ij) * p_d_ij) ^ (-sigma) / P_d_j ^ (1 - sigma);
c_ii = (w_i * L_i + T_i) * ((1 + t_d_ii) * p_d_ii) ^ (-sigma) / P_d_i ^ (1 - sigma);


%% Calculate and compare summary statistics

Omega_ii = M_u_i * M_d_i * p_u_ii * x_ii / (M_d_i * (p_d_ij * c_ij + p_d_ii * c_ii));
Omega_jj = M_u_j * M_d_j * p_u_jj * x_jj / (M_d_j * (p_d_ji * c_ji + p_d_jj * c_jj));
Omega_ij = M_u_i * M_d_j * p_u_ij * x_ij / (M_d_j * (p_d_ji * c_ji + p_d_jj * c_jj));
Omega_ji = M_u_j * M_d_i * p_u_ji * x_ji / (M_d_i * (p_d_ij * c_ij + p_d_ii * c_ii));

b_i_i = M_d_i * p_d_ii * c_ii / (w_i * L_i);
b_i_j = M_d_j * p_d_ji * c_ji / (w_i * L_i);

lambda_d_i = M_d_i * (p_d_ii * c_ii + p_d_ij * c_ij) / (w_i * L_i);

%% Calculate welfare
T_ij=cal_T_ij(t_d_ij,t_u_ij,v_d_ji,v_u_ji,...
                        M_d_i,M_d_j,M_u_i,M_u_j,...
                        x_ij,x_ji,c_ij,c_ji,...
                        p_d_ij,p_d_ji,p_u_ij,p_u_ji);
T_ii=cal_T_ij(t_d_ii,t_u_ii,v_d_ii,v_u_ii,...
                        M_d_i,M_d_i,M_u_i,M_u_i,...
                        x_ii,x_ii,c_ii,c_ii,...
                        p_d_ii,p_d_ii,p_u_ii,p_u_ii);
T_ji=cal_T_ij(t_d_ji,t_u_ji,v_d_ij,v_u_ij,...
                        M_d_j,M_d_i,M_u_j,M_u_i,...
                        x_ji,x_ij,c_ji,c_ij,...
                        p_d_ji,p_d_ij,p_u_ji,p_u_ij);
T_jj=cal_T_ij(t_d_jj,t_u_jj,v_d_jj,v_u_jj,...
                        M_d_j,M_d_j,M_u_j,M_u_j,...
                        x_jj,x_jj,c_jj,c_jj,...
                        p_d_jj,p_d_jj,p_u_jj,p_u_jj);


U_i = cal_U_i(w_i,L_i,T_ii,T_ji,P_d_i);
U_j = cal_U_i(w_j,L_j,T_jj,T_ij,P_d_j);

stats = [Omega_ii, Omega_ji, Omega_jj, Omega_ij, ...
        b_i_i, b_i_j, ...
        lambda_d_i,...
        U_i, U_j];


% Calculate tax revenues
function T_ij=cal_T_ij(t_d_ij,t_u_ij,v_d_ji,v_u_ji,...
                        M_d_i,M_d_j,M_u_i,M_u_j,...
                        x_ij,x_ji,c_ij,c_ji,...
                        p_d_ij,p_d_ji,p_u_ij,p_u_ji)
    T_ij=t_d_ij*M_d_i*c_ij*p_d_ij+...
        t_u_ij*M_d_j*M_u_i*x_ij*p_u_ij-...
        v_d_ji*M_d_j*c_ji*p_d_ji-...
        v_u_ji*M_d_i*M_u_j*x_ji*p_u_ji;
end
% Calculate household utility
function U_i=cal_U_i(w_i,L_i,T_ii,T_ji,P_d_i)
    U_i=(w_i*L_i+T_ii+T_ji)/P_d_i;
end