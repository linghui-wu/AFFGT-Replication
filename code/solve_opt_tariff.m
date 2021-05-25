%% Solve the optimal tarrifs
function opt_u = solve_opt_tariff(x,params)
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
alpha_u = params.alpha_u;
alpha_d = params.alpha_d;

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
t_d_ii = 0; t_d_ij = 0;  % Local tax on downstream goods in country i
t_d_jj = 0;  % Local tax on downstream goods in country j
t_u_ii = 0; t_u_ij = 0;  % Local tax on upstream goods in country i
t_u_jj = 0;  % Local tax on upstream goods in country j

% Subsidies (v_s_ji denotes subsidies to goods from country j sector s to country i)
v_d_ii = 0; v_d_ij = 0; v_d_jj = 0; v_d_ji = 0;
v_u_ii = 0; v_u_ij = 0; v_u_jj = 0; v_u_ji = 0;


%% Equilibrium equations
% Marginal cost of upstream sector in country i and j
mc_u_i=cal_mc_u_i(alpha_u,A_u_i,w_i);
mc_u_j=cal_mc_u_i(alpha_u,A_u_j,x(1)^2);

% Price in upstream sectors
p_u_ij=cal_p_s_ij(mu_u,tau_u,mc_u_i,v_u_ij);
p_u_ji=cal_p_s_ij(mu_u,tau_u,mc_u_j,v_u_ji);
p_u_ii=cal_p_s_ij(mu_u,1,mc_u_i,v_u_ii);
p_u_jj=cal_p_s_ij(mu_u,1,mc_u_j,v_u_jj);

% Price index in upstream sectors
P_u_ij=cal_P_s_ij(x(2)^2,theta,t_u_ij,p_u_ij);
P_u_ii=cal_P_s_ij(x(2)^2,theta,t_u_ii,p_u_ii);
P_u_ji=cal_P_s_ij(x(3)^2,theta,x(8)^2,p_u_ji);
P_u_jj=cal_P_s_ij(x(3)^2,theta,t_u_jj,p_u_jj);

% Price index in upstream sector for country i and country j
P_u_i=cal_P_s_i(P_u_ii,theta,P_u_ji);
P_u_j=cal_P_s_i(P_u_jj,theta,P_u_ij);

% Marginal cost of downstream sectors in country i and j
mc_d_i=cal_mc_d_i(alpha_d,A_d_i,w_i,P_u_i);
mc_d_j=cal_mc_d_i(alpha_d,A_d_j,x(1)^2,P_u_j);

% Price in downstream sectors
p_d_ij=cal_p_s_ij(mu_d,tau_d,mc_d_i,v_d_ij);
p_d_ji=cal_p_s_ij(mu_d,tau_d,mc_d_j,v_d_ji);
p_d_ii=cal_p_s_ij(mu_d,1,mc_d_i,v_d_ii);
p_d_jj=cal_p_s_ij(mu_d,1,mc_d_j,v_d_jj);

% Price index in downstream sectors
P_d_ij=cal_P_s_ij(x(4)^2,sigma,t_d_ij,p_d_ij);
P_d_ji=cal_P_s_ij(x(5)^2,sigma,x(9)^2,p_d_ji);
P_d_ii=cal_P_s_ij(x(4)^2,sigma,t_d_ii,p_d_ii);
P_d_jj=cal_P_s_ij(x(5)^2,sigma,t_d_jj,p_d_jj);

% Price index in downstream sectors in country i and j
P_d_i=cal_P_s_i(P_d_ii,sigma,P_d_ji);
P_d_j=cal_P_s_i(P_d_jj,sigma,P_d_ij);

% Free entry
y_d_i=cal_y_s_i(sigma,f_d); y_d_j=y_d_i;
y_u_i=cal_y_s_i(theta,f_u); y_u_j=y_u_i;

% Quantities in upstream sectors
Q_u_ij=cal_Q_u_ij(alpha_d,mc_d_j,f_d,y_d_j,P_u_j,P_u_ij,theta);
Q_u_ji=cal_Q_u_ij(alpha_d,mc_d_i,f_d,y_d_i,P_u_i,P_u_ji,theta);
Q_u_ii=cal_Q_u_ij(alpha_d,mc_d_i,f_d,y_d_i,P_u_i,P_u_ii,theta);
Q_u_jj=cal_Q_u_ij(alpha_d,mc_d_j,f_d,y_d_j,P_u_j,P_u_jj,theta);

% Output in upstream sectors
x_ij=cal_x_ij(Q_u_ij,t_u_ij,p_u_ij,P_u_ij,theta);
x_ji=cal_x_ij(Q_u_ji,x(8)^2,p_u_ji,P_u_ji,theta);
x_ii=cal_x_ij(Q_u_ii,t_u_ii,p_u_ii,P_u_ii,theta);
x_jj=cal_x_ij(Q_u_jj,t_u_jj,p_u_jj,P_u_jj,theta);

% Output in downstream sectors, currently no tax revenues included
c_ij=cal_c_ij(x(1)^2,L_j,x(7)^2,t_d_ij,p_d_ij,sigma,P_d_j);
c_ji=cal_c_ij(w_i,L_i,x(6)^2,x(9)^2,p_d_ji,sigma,P_d_i);
c_ii=cal_c_ij(w_i,L_i,x(6)^2,t_d_ii,p_d_ii,sigma,P_d_i);
c_jj=cal_c_ij(x(1)^2,L_j,x(7)^2,t_d_jj,p_d_jj,sigma,P_d_j);

%% Equilibrium constraints
% Budget balance
T_ii=cal_T_ij(t_d_ii,t_u_ii,v_d_ii,v_u_ii,x(4)^2,x(4)^2,x(2)^2,x(2)^2,...
            x_ii,x_ii,c_ii,c_ii,p_d_ii,p_d_ii,p_u_ii,p_u_ii);
T_ji=cal_T_ij(x(9)^2,x(8)^2,v_d_ij,v_u_ij,x(5)^2,x(4)^2,x(3)^2,x(2)^2,...
            x_ji,x_ij,c_ji,c_ij,p_d_ji,p_d_ij,p_u_ji,p_u_ij);

%% Function output
opt_u=cal_U_i(w_i,L_i,T_ii,T_ji,P_d_i);

end


%% Define auxilliary functions

% Calculate marginal costs
function mc_u_i=cal_mc_u_i(alpha_u,A_u_i,w_i)
    alpha_u_bar=1/(alpha_u^alpha_u*(1-alpha_u)^(1-alpha_u));
    mc_u_i=alpha_u_bar/A_u_i*w_i^alpha_u;
end

function mc_d_i=cal_mc_d_i(alpha_d,A_d_i,w_i,P_u_i)
    alpha_d_bar=1/(alpha_d^alpha_d*(1-alpha_d)^(1-alpha_d));
    mc_d_i=alpha_d_bar/A_d_i*w_i^alpha_d*P_u_i^(1-alpha_d);
end

% Calculate price in sector s
function p_s_ij=cal_p_s_ij(mu_s,tau_s,mc_s_i,v_s_ij)
    p_s_ij=mu_s*tau_s*mc_s_i/(1+v_s_ij);
end

% Calculate price index in sector s
function P_s_ij=cal_P_s_ij(M_s_i,e,t_s_ij,p_s_ij)
    P_s_ij=M_s_i^(1/(1-e))*(1+t_s_ij)*p_s_ij;
end

function P_s_i=cal_P_s_i(P_s_ii,e,P_s_ji)
    P_s_i=(P_s_ii^(1-e)+P_s_ji^(1-e))^(1/(1-e));
end

% Calculate quantities in the upstream sector
function Q_u_ij=cal_Q_u_ij(alpha_d,mc_d_j,f_d,y_d_j,P_u_j,P_u_ij,theta)
    Q_u_ij=(1-alpha_d)*mc_d_j*(f_d+y_d_j)/P_u_j*(P_u_ij/P_u_j)^-theta;
end

% Calculate output in the upstream sector
function x_ij=cal_x_ij(Q_u_ij,t_u_ij,p_u_ij,P_u_ij,theta)
    x_ij=Q_u_ij*((1+t_u_ij)*p_u_ij/P_u_ij)^-theta;
end

% Calculate output in the downstream sector
function c_ij=cal_c_ij(w_j,L_j,T_j,t_d_ij,p_d_ij,sigma,P_d_j)
    c_ij=(w_j*L_j+T_j)*((1+t_d_ij)*p_d_ij)^-sigma/P_d_j^(1-sigma);
end

% Define free-entry condition
function y_s_i=cal_y_s_i(e,f_i)
    y_s_i=(e-1)*f_i;
end

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
