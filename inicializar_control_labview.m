% INICIALIZACIÓN DEL BUCLE DE CONTROL - IMPLEMENTACIÓN EN LABVIEW
%
% Este script ejecuta la inicialización del observador y el controlador
% implementados en el sistema LabVIEW

import casadi.*

%% ========================================================================
%% CONFIGURACIÓN DEL SISTEMA
%% ========================================================================

nx      = 6;    % Nº estados
nu      = 3;    % Nº entradas
dt_ctrl = 5;    % Periodo control [s]

%% ========================================================================
%% DEFINICIÓN DE PARÁMETROS DEL MODELO DEL CONTROLADOR
%% ========================================================================

p_ctrl.r_D1     = 0.108;    % Radio tanque D1 (m)
p_ctrl.h_D1     = 0.5;      % Altura tanque D1 (m) 
p_ctrl.r_C1     = 0.08;     % Radio tanque C1 (m)
p_ctrl.h_C1     = 0.5;      % Altura tanque C1 (m)
p_ctrl.rho_w    = 1000;     % Densidad agua (kg/m3)
p_ctrl.Cp_w     = 4.186;    % Calor específico agua (kJ/(ºC*kg))
p_ctrl.T_amb    = T_amb;    % Temperatura ambiente (ºC)
p_ctrl.M_w_refr = 3;        % Masa agua disipador (kg)
p_ctrl.M_w_tot  = 13.7;     % Masa agua disipador (kg)
p_ctrl.W_resist = 1;        % Potencia resistencia (kW)

param_modelo_ctrl = struct2array(p_ctrl)';

%% ========================================================================
%% DEFINICIÓN DE CONSTANTES Y PARÁMETROS DE SIMULACIÓN
%% ========================================================================

% Parámetros geométricos y físicos del sistema
p_r_D1    = param_modelo_ctrl(1);   % Radio tanque D1 [m]
p_h_D1    = param_modelo_ctrl(2);   % Altura tanque D1 [m] 
p_r_C1    = param_modelo_ctrl(3);   % Radio tanque C1 [m]
p_h_C1    = param_modelo_ctrl(4);   % Altura tanque C1 [m]
p_rho_w   = param_modelo_ctrl(5);   % Densidad agua [kg/m³]
p_T_amb   = param_modelo_ctrl(7);   % Temperatura ambiente [°C]
p_M_w_tot = param_modelo_ctrl(9);   % Masa agua total D1 + C1 [kg]

% Áreas transversales de tanques
p_A_C1 = pi*p_r_C1^2;   % [m²]
p_A_D1 = pi*p_r_D1^2;   % [m²]

% Parámetros de simulación
t_tot = 3600;                       % Total de tiempo de simulación [s]
N_m   = floor((t_tot + 1)/dt_ctrl); % Número de iteraciones [p]

%% ========================================================================
%% CONFIGURACIÓN DEL ESTIMADOR MHE
%% ========================================================================

% Matriz de pesos para ajuste a mediciones (normalizado)
W_meas = diag([0.1, 0.1, 1, 1, 1]./[15, 100, p_h_C1, 80, 80]);
%             [Q_1 a_v_2 h_w_C1 T_in_D1 T_out_D1]

% Matriz de pesos para regularización de estados (normalizado)
Q_last = diag([2, 2, 2, 2, 2, 2]./[p_h_C1, 0.5, 1, 80, 80, 80]);
%             [h_w_C1 UA d_UA T_in_D1 T_out_D1 T_out_C1]

%% ========================================================================
%% CONFIGURACIÓN DEL CONTROLADOR MPC
%% ========================================================================

% Matriz de pesos por error de seguimiento (h_w_D1 y T_out_D1) (normalizado)
% W_track = diag([1, 10]./[p_h_D1, 80]);
%                [h_w_D1 T_out_D1]


% Matriz de pesos por variación de control (normalizado)
% R_ctrl  = 0.01;                       % Peso global vs W_track
W_delta = diag([3, 3, 0.5]./[1, 1, 1]); % Pesos individuales
%              [a_v_1 a_v_2 w_vent]

%% ========================================================================
%% CONSTRUCCIÓN DE MODELO Y SOLVER DEL OBSERVADOR
%% ========================================================================

% Crear modelo dinámico simbólico
ode_mhe = modelo_mhe;

% Construir integrador numérico
simulador_mhe = constructor_integrador(ode_mhe, dt_ctrl);

% Construir solver del observador MHE con Multiple Shooting
[solver_mhe, v_est_min, v_est_max, g_est_min, g_est_max] = ...
    constructor_mhe_MS(N_e, W_meas, Q_last, param_modelo_ctrl, simulador_mhe);

%% ========================================================================
%% CONSTRUCCIÓN DE MODELO Y SOLVER DEL CONTROLADOR
%% ========================================================================

% Crear modelo dinámico simbólico
ode_mpc = modelo_mpc;

% Construir integrador numérico
simulador_mpc = constructor_integrador(ode_mpc, dt_ctrl);

% Construir solver del controlador MPC con Single Shooting
[solver_mpc, u_ctrl_min, u_ctrl_max, g_ctrl_min, g_ctrl_max] = ...
    constructor_mpc_SS(N_h, N_c, W_delta, param_modelo_ctrl, simulador_mpc);

%% ========================================================================
%% INICIALIZACIÓN DE OBSERVADOR Y CONTROLADOR
%% ========================================================================

% Estado inicial
x_k = [y_meas_k(1);     % h_w_C1 desde medición [m] 
       0;               % UA inicial [kW/°C]
       0;               % d_UA inicial [kW/(°C·s)]
       y_meas_k(2);     % T_in_D1 desde medición [°C]
       y_meas_k(3);     % T_out_D1 desde medición [°C]
       y_meas_k(3)];    % T_out_C1 inicial [°C]

% Inicializar variables del observador MHE
[v_est_seq, y_meas_buff, u_buff, lam_g_est, lam_x_est] = ...
    inicializar_observador(N_e, size(g_est_min, 1), x_k, y_meas_k, u_meas_k);

% Inicializar variables del controlador MPC
[u_ctrl_seq, x_pred, lam_g_ctrl, lam_x_ctrl] = ...
    inicializar_controlador(N_c, size(g_ctrl_min, 1), x_k, u_k);

% Status de control
status_est_k  = '';
status_ctrl_k = '';