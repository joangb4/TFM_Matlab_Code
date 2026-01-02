% SIMULACIÓN Y VALIDACIÓN DEL MODELO
% 
% Este script ejecuta la simulación completa del modelo del sistema
% usando datos experimentales reales del proceso

clear; clc;
import casadi.*

% ========================================================================
% CARGA Y PREPARACIÓN DE DATOS EXPERIMENTALES
% ========================================================================

fprintf('=== SIMULACIÓN Y VALIDACIÓN DEL MODELO ===\n');
fprintf('Cargando datos experimentales...\n');

% Cargar datos de ensayo experimental
load('Experimento_07_01_2025.mat');

exp_t    = exp(:, 1);       % Datos experimentales de tiempo [s]       
exp_u    = exp(:, 2:4);     % Datos experimentales de acciones de control [a_v_1, a_v_2, w_vent]
exp_niv  = exp(:, 5:9);     % Datos experimentales de balance volumétrico [Q_1_Lmin, h_w_D1, h_w_D1_per, h_w_C1, h_w_C1_per]
exp_temp = exp(:, 10:11);   % Datos experimentales de balance energético [T_in_D1, T_out_D1]

dt_exp = 1;             % Periodo de muestreo del experimento [s/p]
t_tot  = exp_t(end);    % Tiempo total experimento [s]

% Histórico de matrices de entrada y salida medida del sistema
u_hist      = exp_u;                                            % [a_v_1, a_v_2, w_vent]
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

% Periodo de simulación [s/p]
dt_sim = 1;

% Número de iteraciones [p]
N_m   = floor((t_tot + 1)/dt_sim);  

fprintf('Datos cargados: %d muestras, periodo de simulación: %d s\n', N_m, dt_sim);

% ========================================================================
% CONSTRUCCIÓN DE MODELO Y SOLVER
% ========================================================================

fprintf('Construyendo modelo dinámico y solver...\n');

% Crear modelo dinámico simbólico
ode = modelo_mpc;

% Construir simulador
simulador = constructor_integrador(ode, dt_sim);

% ========================================================================
% INICIALIZACIÓN DEL SIMULADOR
% ========================================================================

fprintf('Inicializando simulador con condiciones iniciales...\n');

% Condiciones iniciales del sistema
u_k      = u_hist(1, :)';       % Primera entrada
u_resist = 1;                   % Resistencia calefactora ON
y_meas_k = y_meas_hist(1, :)';  % Primera medición

% Estado inicial estimado
x_k = [y_meas_k(1);     % h_w_C1 desde medición [m] 
       0;               % UA inicial [kW/°C]
       0;               % d_UA inicial [kW(°C·s)]
       y_meas_k(2);     % T_in_D1 desde medición [°C]
       y_meas_k(3);     % T_out_D1 desde medición [°C]
       y_meas_k(3)];    % T_out_C1 inicial [°C]

fprintf(['Estado inicial: h_C1 = %.3f m, UA = %.3f kW/°C, d_UA = %.3f kW/(°C·s)\n' ...
         '                T_in_D1 = %.1f °C, T_out_D1 = %.1f °C, T_out_C1 = %.1f °C\n'], ...
        x_k(1), x_k(2), x_k(3), x_k(4), x_k(5), x_k(6));

% ========================================================================
% BUCLE PRINCIPAL DE SIMULACIÓN
% ========================================================================

fprintf('\n=== INICIANDO SIMULADOR ===\n');

% Preparación de variables de resultados
sim_u           = zeros(N_m, nu);   % Acciones de control    
sim_x           = zeros(N_m, nx);   % Estados
sim_t           = zeros(N_m, 1);    % Tiempos
sim_tiempo_comp = zeros(N_m, 1);    % Tiempos de cómputo

% Bucle principal
for k = 1:N_m
    
    % Indexado datos experimentales
    idx_sim = (k - 1)*dt_sim + 1;
    
    % Obtener nuevas mediciones y entradas
    u_k = u_hist(idx_sim, :)';
    
    % Ejecutar paso de simulación
    tic;

    x_kp1 = simulador_1p(simulador, param_modelo, x_k, u_k, u_resist);
    
    t_k = toc;

    % Almacenar resultados
    sim_u(k, :)        = u_k';
    sim_x(k, :)        = x_k';
    sim_t(k)           = (k - 1)*dt_sim;
    sim_tiempo_comp(k) = t_k;
    
    % Actualización de estados
    x_k = x_kp1;
    
    % Mostrar progreso
    if mod(k, 20) == 0 || k <= 5
        fprintf('Iter %3d → Tiempo: %6.2f ms\n', ...
                k, t_k*1000);
    end
end

fprintf('\n=== SIMULACIÓN COMPLETADA ===\n');
fprintf('Total iteraciones: %d\n', k);
fprintf('Tiempo promedio: %.2f ms\n', mean(sim_tiempo_comp)*1000);
fprintf('Tiempo máximo: %.2f ms\n', max(sim_tiempo_comp)*1000);
