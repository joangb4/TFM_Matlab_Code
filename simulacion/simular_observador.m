% SIMULACIÓN Y VALIDACIÓN DEL OBSERVADOR
% 
% Este script ejecuta la simulación completa del sistema de observación
% usando datos experimentales reales del proceso

clear; clc;
import casadi.*

% ========================================================================
% CARGA Y PREPARACIÓN DE DATOS EXPERIMENTALES
% ========================================================================

fprintf('=== SIMULACIÓN Y VALIDACIÓN DEL OBSERVADOR ===\n');
fprintf('Cargando datos experimentales...\n');

% Cargar datos de ensayo experimental
load('Experimento_07_01_2025.mat');

exp_t    = exp(:, 1);       % Datos experimentales de tiempo [s]       
exp_u    = exp(:, 2:4);     % Datos experimentales de acciones de control [a_v_1, a_v_2, w_vent]
exp_niv  = exp(:, 5:9);     % Datos experimentales de balance volumétrico [Q_1_Lmin, h_w_D1, h_w_D1_per, h_w_C1, h_w_C1_per]
exp_temp = exp(:, 10:11);   % Datos experimentales de balance energético [T_in_D1, T_out_D1]

dt_exp = 1;             % Periodo de muestreo del experimento [s/p]
t_tot  = exp_t(end);    % Tiempo total experimento [s]

% Matrices de entrada y salida medida del observador
u_hist      = [exp_niv(:, 1), exp_u(:, 2:3)];                   % [Q_1_Lmin, a_v_2, w_vent]
y_meas_hist = [exp_niv(:, 4), exp_temp(:, 1), exp_temp(:, 2)];  % [h_C1, T_in_D1, T_out_D1]

% ========================================================================
% DEFINICIÓN DE PARÁMETROS DEL MODELO FÍSICO
% ========================================================================

fprintf('Configurando parámetros del modelo físico...\n');

% Parámetros geométricos y físicos del sistema
p.r_D1     = 0.108; % Radio tanque D1 [m]
p.h_D1     = 0.5;   % Altura tanque D1 [m] 
p.r_C1     = 0.08;  % Radio tanque C1 [m]
p.h_C1     = 0.5;   % Altura tanque C1 [m]
p.rho_w    = 1000;  % Densidad agua [kg/m³]
p.Cp_w     = 4.186; % Calor específico agua [kJ/(°C·kg)]
p.T_amb    = 20;    % Temperatura ambiente [°C]
p.M_w_refr = 3;     % Masa agua radiador [kg]
p.M_w_tot  = 13.7;  % Masa agua total D1 + C1 [kg]
p.W_resist = 1;     % Potencia resistencia [kW]

param_modelo = struct2array(p)';

% Dimensiones del problema
nx = 6; % Número de estados
nu = 3; % Número de entradas

% ========================================================================
% DEFINICIÓN DE CONSTANTES Y PARÁMETROS DE SIMULACIÓN
% ========================================================================

fprintf('Configurando constantes y parámetros de simulación...\n');

% Parámetros geométricos y físicos del sistema
p_r_D1    = 0.108;  % Radio tanque D1 [m]
p_h_D1    = 0.5;    % Altura tanque D1 [m] 
p_r_C1    = 0.08;   % Radio tanque C1 [m]
p_h_C1    = 0.5;    % Altura tanque C1 [m]
p_rho_w   = 1000;   % Densidad agua [kg/m³]
p_T_amb   = 20;     % Temperatura ambiente [°C]
p_M_w_tot = 13.7;   % Masa agua total D1 + C1 [kg]

% Áreas transversales de tanques
p_A_C1 = pi*p_r_C1^2;   % [m²]
p_A_D1 = pi*p_r_D1^2;   % [m²]

% Parámetros de simulación
dt_sim = 5;                         % Periodo de simulación [s/p]
N_m    = floor((t_tot + 1)/dt_sim); % Número de iteraciones [p]

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
% CONSTRUCCIÓN DE MODELO Y SOLVER
% ========================================================================

fprintf('Construyendo modelo dinámico y solver...\n');

% Crear modelo dinámico simbólico
ode_mhe = modelo_mhe();

% Periodo de simulación del observador
dt_mhe = 5;

% Construir integrador numérico
simulador_mhe = constructor_integrador(ode_mhe, dt_mhe);

% Construir solver del observador MHE con Multiple Shooting
[solver_mhe, v_est_min, v_est_max, g_est_min, g_est_max] = ...
    constructor_mhe_MS(N_e, W_meas, W_reg, param_modelo, simulador_mhe);

fprintf('Solver MHE construido exitosamente\n');

% ========================================================================
% INICIALIZACIÓN DEL ESTIMADOR
% ========================================================================

fprintf('Inicializando estimador con condiciones iniciales...\n');

% Condiciones iniciales del sistema
u_k      = u_hist(1, :)';       % Primera entrada
u_resist = 1;                   % Resistencia calefactora ON
y_meas_k = y_meas_hist(1, :)';  % Primera medición

% Estado inicial estimado
x_k = [y_meas_k(1);     % h_w_C1 desde medición [m] 
       0;               % UA inicial [kW/°C]
       0;               % d_UA inicial [kW/(°C·s)]
       y_meas_k(2);     % T_in_D1 desde medición [°C]
       y_meas_k(3);     % T_out_D1 desde medición [°C]
       y_meas_k(3)];    % T_out_C1 inicial [°C]

% Inicializar variables del estimador
[v_est_seq, y_meas_buff, u_buff, lam_g_mhe, lam_x_mhe] = ...
    inicializar_observador(N_e, size(g_est_min, 1), x_k, y_meas_k, u_k);

fprintf(['Estado inicial: h_C1 = %.3f m, UA = %.3f kW/°C, d_UA = %.3f kW/(°C·s)\n' ...
         '                T_in_D1 = %.1f °C, T_out_D1 = %.1f °C, T_out_C1 = %.1f °C\n'], ...
        x_k(1), x_k(2), x_k(3), x_k(4), x_k(5), x_k(6));

% ========================================================================
% BUCLE PRINCIPAL DE ESTIMACIÓN
% ========================================================================

fprintf('\n=== INICIANDO ESTIMACIÓN MHE ===\n');

% Preparación de variables de resultados
est_u           = zeros(N_m, 2);    % Acciones de control
est_x           = zeros(N_m, nx);   % Estados
est_t           = zeros(N_m, 1);    % Tiempos
est_tiempo_comp = zeros(N_m, 1);    % Tiempos de cómputo
est_status      = cell(N_m, 1);     % Estado del estimador

% Bucle principal
for k = 1:N_m

    % Indexado datos experimentales
    idx_sim = (k - 1)*dt_sim + 1;
    
    % Obtener nuevas mediciones y entradas
    u_k      = u_hist(idx_sim, :)';
    y_meas_k = y_meas_hist(idx_sim, :)';
    
    % Ejecutar paso de estimación MHE
    tic;

    [v_est_seq, x_est_k, u_est_k, lam_g_mhe, lam_x_mhe, ...
     y_meas_buff, u_buff, status_mhe_k, n_iter_mhe_k] = ...
        estimar_estado(N_e, solver_mhe, v_est_seq, lam_g_mhe, lam_x_mhe, ...
                      v_est_min, v_est_max, g_est_min, g_est_max, ...
                      y_meas_buff, y_meas_k, u_buff, u_k, u_resist);
    t_k = toc;

    % Almacenar resultados
    est_u(k, :)        = u_est_k';
    est_x(k, :)        = x_est_k';    
    est_t(k)           = (k - 1)*dt_mhe;
    est_tiempo_comp(k) = t_k;
    est_status{k}      = status_mhe_k;
    
    % Mostrar progreso
    if mod(k, 20) == 0 || k <= 5
        fprintf('Iter %3d → MHE IPOPT iters: %3d, tiempo: %6.2f ms, status: %s\n', ...
                k, n_iter_mhe_k, t_k*1000, status_mhe_k);
    end
end

fprintf('\n=== ESTIMACIÓN COMPLETADA ===\n');
fprintf('Total iteraciones MHE: %d\n', k);
fprintf('Tiempo de iteración promedio: %.2f ms\n', mean(est_tiempo_comp)*1000);
fprintf('Tiempo de iteración máximo: %.2f ms\n', max(est_tiempo_comp)*1000);

% Análisis de convergencia
status_exitosos = sum(contains(est_status, 'Solve_Succeeded'));
fprintf('Convergencia exitosa: %.1f%% (%d/%d)\n', ...
        status_exitosos/k*100, status_exitosos, k);
