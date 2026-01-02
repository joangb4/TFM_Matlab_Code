import casadi.*
cd 'D:\MATLAB\TFM_data'

% ========================================================================
% ENTRADAS DEL BUCLE DE CONTROL
% ========================================================================

% Además de N_h, N_c, N_e, T_amb y dt_ctrl
u_k      = u_k(:);      % Acción de control
ref_k    = ref_k(:);    % Referencias horizonte
y_meas_k = y_meas_k(:); % Salidas medidas
u_meas_k = u_meas_k(:); % Acción de control medida para observador

% Ensamble de matriz de pesos por error de seguimiento
w_track = [w_track_h, w_track_T];

if k == 1
    addpath('D:\MATLAB\casadi-3.6.7-windows64-matlab2018b')
    % Inicializar sistema
    inicializar_control_labview;
end

obtener_control_labview;

% ========================================================================
% SALIDAS DEL BUCLE DE CONTROL A LABVIEW
% ========================================================================

u_k = u_k';
u_est_k = u_est_k';
x_est_k = x_est_k';
