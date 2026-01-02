% SIMULACIÓN Y VALIDACIÓN DEL CONTROLADOR
%
% Este script ejecuta la simulación completa del sistema de control
% a partir de un estado inicial y unas consignas impuestas

clear; clc;
import casadi.*

fprintf('=== SIMULACIÓN Y VALIDACIÓN DEL CONTROLADOR ===\n');

% ========================================================================
% CONFIGURACIÓN DEL SISTEMA
% ========================================================================

fprintf('Configurando el sistema...\n');

nx      = 6;    % Nº estados
nu      = 3;    % Nº entradas
dt_ctrl = 5;    % Periodo control [s]

% ========================================================================
% DEFINICIÓN DE PARÁMETROS DEL MODELO FÍSICO
% ========================================================================

fprintf('Configurando parámetros del modelo físico...\n');

% Parámetros geométricos y físicos del sistema
p.r_D1     = 1.2*0.108; % Radio tanque D1 [m]
p.h_D1     = 1.2*0.5;   % Altura tanque D1 [m] 
p.r_C1     = 0.8*0.08;  % Radio tanque C1 [m]
p.h_C1     = 0.8*0.5;   % Altura tanque C1 [m]
p.rho_w    = 1000;      % Densidad agua [kg/m³]
p.Cp_w     = 4.186;     % Calor específico agua [kJ/(°C·kg)]
p.T_amb    = 24;        % Temperatura ambiente [°C]
p.M_w_refr = 0.8*3;     % Masa agua radiador [kg]
p.M_w_tot  = 1.2*13.7;  % Masa agua total D1 + C1 [kg]
p.W_resist = 0.8*1;     % Potencia resistencia [kW]

param_modelo = struct2array(p)';

% ========================================================================
% DEFINICIÓN DE PARÁMETROS DEL MODELO DEL CONTROLADOR
% ========================================================================

fprintf('Configurando parámetros del modelo del controlador...\n');

p_ctrl.r_D1     = 0.108;    % Radio tanque D1 (m)
p_ctrl.h_D1     = 0.5;      % Altura tanque D1 (m) 
p_ctrl.r_C1     = 0.08;     % Radio tanque C1 (m)
p_ctrl.h_C1     = 0.5;      % Altura tanque C1 (m)
p_ctrl.rho_w    = 1000;     % Densidad agua (kg/m3)
p_ctrl.Cp_w     = 4.186;    % Calor específico agua (kJ/(ºC*kg))
p_ctrl.T_amb    = 23;       % Temperatura ambiente (ºC)
p_ctrl.M_w_refr = 3;        % Masa agua disipador (kg)
p_ctrl.M_w_tot  = 13.7;     % Masa agua disipador (kg)
p_ctrl.W_resist = 1;        % Potencia resistencia (kW)

param_modelo_ctrl = struct2array(p_ctrl)';

% ========================================================================
% DEFINICIÓN DE CONSTANTES Y PARÁMETROS DE SIMULACIÓN
% ========================================================================

fprintf('Configurando constantes y parámetros de simulación...\n');

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

% ========================================================================
% CONFIGURACIÓN DEL ESTIMADOR MHE
% ========================================================================

fprintf('Configurando estimador MHE...\n');

% Parámetros del estimador
N_e = 6;    % Horizonte de estimación [p]

% Matriz de pesos para ajuste a mediciones (normalizado)
W_meas = diag([0.1, 0.1, 1, 1, 1]./[15, 100, p_h_C1, 80, 80]);
%             [Q_1 a_v_2 h_w_C1 T_in_D1 T_out_D1]

% Matriz de pesos para regularización de estados (normalizado)
W_reg = diag([2, 2, 2, 2, 2, 2]./[p_h_C1, 0.5, 1, 80, 80, 80]);
%            [h_w_C1 UA d_UA T_in_D1 T_out_D1 T_out_C1]

fprintf('Horizonte MHE: %d pasos, Pesos configurados\n', N_e);

% ========================================================================
% CONFIGURACIÓN DEL CONTROLADOR MPC
% ========================================================================

fprintf('Configurando controlador MPC...\n');

% Parámetros del controlador
N_h = 20;   % Horizonte de predicción [p]
N_c = 5;    % Horizonte de control [p]

% Matriz de pesos por error de seguimiento (h_w_D1 y T_out_D1) (normalizado)
w_track = [1, 10];
%         [h_w_D1 T_out_D1]

% Matriz de pesos por variación de control (normalizado)
R_ctrl  = 0.02;                         % Peso global vs W_track
W_delta = diag([3, 3, 0.5]./[1, 1, 1]); % Pesos individuales
%              [a_v_1 a_v_2 w_vent]

fprintf('Horizonte MPC: %d pasos de predicción, %d pasos de control, Pesos configurados\n', N_h, N_c);

% ========================================================================
% CONSTRUCCIÓN DE MODELO Y SOLVER DEL OBSERVADOR
% ========================================================================

fprintf('Construyendo modelo dinámico y solver MHE...\n');

% Crear modelo dinámico simbólico
ode_mhe = modelo_mhe;

% Construir integrador numérico
simulador_mhe = constructor_integrador(ode_mhe, dt_ctrl);

% Construir solver del observador MHE con Multiple Shooting
[solver_mhe, v_est_min, v_est_max, g_est_min, g_est_max] = ...
    constructor_mhe_MS(N_e, W_meas, W_reg, param_modelo_ctrl, simulador_mhe);

fprintf('Solver MHE construido exitosamente\n');

% ========================================================================
% CONSTRUCCIÓN DE MODELO Y SOLVER DEL CONTROLADOR
% ========================================================================

fprintf('Construyendo modelo dinámico y solver MPC...\n');

% Crear modelo dinámico simbólico
ode_mpc = modelo_mpc;

% Construir integrador numérico
simulador_mpc = constructor_integrador(ode_mpc, dt_ctrl);

% Construir solver del controlador MPC con Single Shooting
[solver_mpc, u_ctrl_min, u_ctrl_max, g_ctrl_min, g_ctrl_max] = ...
    constructor_mpc_SS(N_h, N_c, W_delta, param_modelo_ctrl, simulador_mpc);

fprintf('Solver MPC construido exitosamente\n');

% ========================================================================
% INICIALIZACIÓN DE OBSERVADOR Y CONTROLADOR
% ========================================================================

fprintf('Inicializando controlador con condiciones iniciales...\n');

% Condiciones iniciales del sistema
u_k            = [50; 42; 17];                      % Primera entrada
u_resist       = 1;                                 % Resistencia calefactora ON
Q_1_Lmin_sim_k = 0.00195*u_k(1)^2 - 0.023*u_k(1);   % [L/min]
y_meas_k       = [0.4*p_h_C1; 25; 25];              % Primera salida medida
u_meas_k       = [Q_1_Lmin_sim_k; u_k(2:3)];        % Primera entrada medida

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

fprintf(['Estado inicial: h_C1 = %.3f m, UA = %.3f kW/°C, d_UA = %.3f kW/(°C·s)\n' ...
         '                T_in_D1 = %.1f °C, T_out_D1 = %.1f °C, T_out_C1 = %.1f °C\n'], ...
        x_k(1), x_k(2), x_k(3), x_k(4), x_k(5), x_k(6));

% ========================================================================
% GENERACIÓN DE REFERENCIAS
% ========================================================================

% Referencias de h_w_D1 en 50% (evita vaciado-llenado como sol)
ref_h_D1 = 50*p_h_D1/100;
ref_h_D1 = repmat(ref_h_D1, N_m + N_h, 1);

% Referencias de T_out_D1
% En escalones
ref_T_out_D1 = vertcat(30*ones(150, 1), ...
                       35*ones(150, 1), ...
                       45*ones(270, 1), ...
                       40*ones(N_m + N_h - 570, 1));

% Simulando proceso
% ref_T_out_D1 = vertcat((30*ones(1, 50))', ...
%                        linspace(30, 47, 300)', ...
%                        (50*ones(1, 150))', ...
%                        linspace(50, 30, 150)', ...
%                        (30*ones(1, max(1, N_m + N_h - 900)))');

% Ensamble de vector referencias
ref          = [ref_h_D1, ref_T_out_D1];
ref_reshaped = reshape(ref.', [], 1);

% ========================================================================
% BUCLE PRINCIPAL DE CONTROL
% ========================================================================

fprintf('\n=== INICIANDO CONTROL MPC ===\n');

% Preparación de variables de resultados
% Variables de estimador
est_u      = zeros(N_m, 2);     % Acciones de control
est_x      = zeros(N_m, nx);    % Estados
est_status = cell(N_m, 1);      % Estado del estimador

% Variables de controlador
ctrl_u           = zeros(N_m, 3);   % Acciones de control  
ctrl_tiempo_comp = zeros(N_m, 1);   % Tiempos de cómputo
ctrl_status      = cell(N_m, 1);    % Estado del controlador
ctrl_ures        = zeros(N_m, 1);   %

% Variables de simulación
sim_Q_1 = zeros(N_m, 1);
sim_x = zeros(N_m, nx); % Estados
sim_t = zeros(N_m, 1);  % Tiempos

% Bucle principal
for k = 1:N_m

    % Estados sensados [h_w_C1, T_in_D1, T_out_D1]
    y_meas_k = [x_k(1); x_k(4); x_k(5)];

    % Acción de control (formato estimador) [Q_1_Lmin, a_v_2, w_vent]
    u_meas_k = [Q_1_Lmin_sim_k; u_k(2:3)];

    % Referencias en k
    idx_ref_k = 2*(k - 1) + 1 : 2*(k - 1) + 2*N_h;
    ref_k     = ref_reshaped(idx_ref_k);

    tic

    [v_est_seq, x_est_k, u_est_k, lam_g_est, lam_x_est, ...
     y_meas_buff, u_buff, status_est_k, n_iter_est_k] = estimar_estado(N_e, solver_mhe, v_est_seq, ...
                                                                       lam_g_est, lam_x_est, ...
                                                                       v_est_min, v_est_max, ...
                                                                       g_est_min, g_est_max, ...
                                                                       y_meas_buff, y_meas_k, ...
                                                                       u_buff, u_meas_k, u_resist);

    [u_ctrl_seq, u_ctrl_k, x_pred, ...
     lam_g_ctrl, lam_x_ctrl, status_ctrl_k, n_iter_ctrl_k] = optimizar_control(solver_mpc, u_ctrl_seq, ...
                                                                               lam_g_ctrl, lam_x_ctrl, ...
                                                                               u_ctrl_min, u_ctrl_max, ...
                                                                               g_ctrl_min, g_ctrl_max, ...
                                                                               x_est_k, u_k, u_resist, ...
                                                                               ref_k, x_pred, y_meas_k, ...
                                                                               simulador_mpc, w_track, ...
                                                                               R_ctrl, param_modelo_ctrl);
    
    t_k = toc;
    
    % Acción de control obtenida
    u_k = u_ctrl_k;

    % Cálculo de caudal 1 (bomba G1)
    Q_1_Lmin_sim_k = 0.00195*u_k(1)^2 - 0.023*u_k(1);   % [L/min]

    % Simulación del comportamiento de la planta
    x_kp1 = simulador_1p(simulador_mpc, param_modelo, x_k, u_k, u_resist);

    % Simulación de ruido blanco en las variables sensadas (simulación de sensores)
    sigma_Q_1Lmin  = 0.15;
    sigma_h_C1     = 0.2/100*p_h_C1;
    sigma_T_in_D1  = 0.07;
    sigma_T_out_D1 = 0.07;

    sigmas = [sigma_Q_1Lmin; sigma_h_C1; sigma_T_in_D1; sigma_T_out_D1];

    meas_k = aplicar_ruido([Q_1_Lmin_sim_k; x_kp1(1); x_kp1(4); x_kp1(5)], sigmas);

    Q_1_Lmin_sim_k = meas_k(1);
    x_kp1(1)       = meas_k(2);
    x_kp1(4)       = meas_k(3);
    x_kp1(5)       = meas_k(4);

    % Almacenar resultados
    est_u(k, :)   = u_est_k';
    est_x(k, :)   = x_est_k';
    est_status{k} = status_est_k;

    ctrl_u(k, :)        = u_k';
    ctrl_tiempo_comp(k) = t_k;
    ctrl_status{k}      = status_ctrl_k;
    ctrl_ures(k)        = u_resist;

    sim_Q_1(k, :) = Q_1_Lmin_sim_k;
    sim_x(k, :)   = x_k';
    sim_t(k)      = (k - 1)*dt_ctrl;

    % Actualización de estados
    x_k = x_kp1;
    
    % Mostrar progreso
    if mod(k, 20) == 0 || k <= 5
        fprintf('Iter %3d → MHE iters: %3d, MHE status: %s, MPC iters: %3d, MPC status: %s, tiempo total: %.2f ms\n', ...
                k, n_iter_est_k, status_est_k, n_iter_ctrl_k, status_ctrl_k, t_k*1000);
    end    
end

fprintf('\n=== CONTROL COMPLETADO ===\n');
fprintf('Total iteraciones: %d\n', k);
fprintf('Tiempo de iteración promedio: %.2f ms\n', mean(ctrl_tiempo_comp)*1000);
fprintf('Tiempo de iteración máximo: %.2f ms\n', max(ctrl_tiempo_comp)*1000);

% Análisis de convergencia (MHE y MPC)
status_exitosos_est = sum(contains(est_status, 'Solve_Succeeded'));
fprintf('Convergencia exitosa (MHE): %.1f%% (%d/%d)\n', ...
        status_exitosos_est/k*100, status_exitosos_est, k);

status_exitosos_ctrl = sum(contains(ctrl_status, 'Solve_Succeeded'));
fprintf('Convergencia exitosa (MPC): %.1f%% (%d/%d)\n', ...
        status_exitosos_ctrl/k*100, status_exitosos_ctrl, k);
