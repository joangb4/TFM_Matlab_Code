function [solver, v_min, v_max, g_min, g_max] = constructor_mhe_MS(N_e, W_meas, W_reg, param_modelo, simulador)
    % CONSTRUCTOR_MHE_MS - Construye optimizador MHE con Multiple Shooting
    %
    % Descripción:
    %   Implementa la formulación de Moving Horizon Estimation (MHE) usando
    %   la técnica de Multiple Shooting para estimación simultánea de estados
    %   y acciones de control inciertas (Q_1_Lmin, a_v_2)
    %
    % Entradas:
    %   N_e          - Horizonte de estimación [p]
    %   W_meas       - Matriz de pesos de mediciones
    %   W_last       - Matriz de pesos de regularización
    %   param_modelo - Parámetros del modelo físico
    %   simulador    - Integrador numérico del modelo
    %
    % Salidas:
    %   solver       - Optimizador IPOPT configurado
    %   v_min        - Límites inferiores de variables a estimar
    %   v_max        - Límites superiores de variables a estimar
    %   g_min        - Límites inferiores de restricciones
    %   g_max        - Límites superiores de restricciones
    
    import casadi.*
    
    % =====================================================================
    % DEFINICIÓN DE PARÁMETROS DEL PROBLEMA
    % =====================================================================
    
    nx = 6; % Número de estados
    nu = 3; % Número de entradas
    
    % Extracción de parámetros físicos del modelo
    p_h_C1  = param_modelo(4);  % Altura tanque C1 [m]
    p_T_amb = param_modelo(7);  % Temperatura ambiente [°C]
    
    % =====================================================================
    % VARIABLES DE DECISIÓN DEL PROBLEMA MHE
    % =====================================================================
    
    % Variables a estimar
    x_est = MX.sym('x_est', nx*(N_e + 1));  % Estados [x0, x1, ..., xN_e]
    u_est = MX.sym('u_est', 2*N_e);         % Entradas inciertas [Q_1_Lmin, a_v_2]
    
    % Parámetros del problema (datos conocidos)
    y_meas_buff = MX.sym('y_meas_buff', 3*N_e);         % Histórico de salidas medidas [h_w_C1, T_in_D1, T_out_D1]
    u_buff      = MX.sym('u_buff', nu*N_e);             % Histórico de entradas [Q_1_Lmin, a_v_2, w_vent]
    u_resist    = MX.sym('u_resist', 1);                % Control resistencia [ON/OFF]
    x_est_last  = MX.sym('x_est_last', nx*(N_e + 1));   % Estimación anterior
    
    % =====================================================================
    % INICIALIZACIÓN DE FUNCIÓN OBJETIVO Y RESTRICCIONES
    % =====================================================================
    
    J = 0;  % Función objetivo
    g = []; % Vector de restricciones
    
    % =====================================================================
    % BUCLE DE HORIZONTE DE ESTIMACIÓN (N_e)
    % =====================================================================
    
    for k = 1:N_e
        
        % Indexación de variables en instante k
        idx_x_k     = (k - 1)*nx + 1 : k*nx;    % Estados x_k
        idx_x_kp1   = k*nx + 1 : (k + 1)*nx;    % Estados x_{k+1}
        idx_u_k     = (k - 1)*nu + 1 : k*nu;    % Entradas u_k
        idx_u_est_k = (k - 1)*2 + 1 : k*2;      % Entradas estimadas
        idx_y_k     = (k - 1)*3 + 1 : k*3;      % Mediciones y_k
        
        % Extracción de variables
        x_k      = x_est(idx_x_k);          % Estado actual
        x_kp1    = x_est(idx_x_kp1);        % Estado siguiente
        u_k      = u_buff(idx_u_k);         % Entrada conocida
        u_est_k  = u_est(idx_u_est_k);      % Entrada estimada [Q_1, a_v_2]
        y_meas_k = y_meas_buff(idx_y_k);    % Medición actual
        
        % Construcción del vector de entrada completo
        u_meas_k     = u_k(1:2);            % Entradas medidas
        u_completo_k = [u_est_k; u_k(3)];   % Entradas estimadas + w_vent [Q_1_est, a_v_2_est, w_vent]
        
        % Simulación del modelo dinámico
        params_sim_k = vertcat(u_completo_k, u_resist, param_modelo);
        sim_res      = simulador('x0', x_k, 'p', params_sim_k);
        x_pred_kp1   = sim_res.xf;  % Estado predicho
        
        % Construcción de salidas estimadas
        x_meas_est_kp1 = [x_kp1(1);    % h_w_C1 estimado
                          x_kp1(4);    % T_in_D1 estimado
                          x_kp1(5)];   % T_out_D1 estimado
        
        % Vectores de mediciones y estimaciones
        v_meas_kp1     = vertcat(u_meas_k, y_meas_k);       % [Q_1_med, a_v_2_med, mediciones]
        v_meas_est_kp1 = vertcat(u_est_k, x_meas_est_kp1);  % [Q_1_est, a_v_2_est, estimaciones]
        
        % Término de ajuste a mediciones en función objetivo
        dv = v_meas_kp1 - v_meas_est_kp1;
        J  = J + dv'*W_meas*dv;
        
        % Restricción de costura dinámica (Multiple Shooting)
        g = [g; x_kp1 - x_pred_kp1];
    end
    
    % =====================================================================
    % TÉRMINO DE REGULARIZACIÓN
    % =====================================================================
    
    % Penalización por cambios de estados respecto a estimación anterior
    dx_last = x_est(1:nx) - x_est_last(nx + 1 : 2*nx);
    J = J + dx_last'*W_reg*dx_last;
    
    % =====================================================================
    % DEFINICIÓN DE LÍMITES Y RESTRICCIONES
    % =====================================================================
    
    % Límites de estados físicos
    x_min = repmat([0.01*p_h_C1;        % h_w_C1 mínimo  
                    -1e-6;              % UA mínimo
                    -inf;               % d_UA mínimo
                    p_T_amb - 5;        % T_in_D1 mínimo
                    p_T_amb - 5;        % T_out_D1 mínimo
                    p_T_amb - 5], ...   % T_out_C1 mínimo
                   (N_e + 1), 1);
    
    x_max = repmat([0.98*p_h_C1;    % h_w_C1 máximo
                    inf;            % UA máximo
                    inf;            % d_UA máximo
                    80;             % T_in_D1 máximo
                    80;             % T_out_D1 máximo
                    80], ...        % T_out_C1 máximo
                   (N_e + 1), 1);
    
    % Límites de entradas estimadas
    u_meas_min = repmat([0;         % Q_1_Lmin mínimo
                         0], ...    % a_v_2 mínimo
                        N_e, 1);
    
    u_meas_max = repmat([15;        % Q_1_Lmin máximo
                         100], ...  % a_v_2 máximo
                        N_e, 1);
    
    % Ensamble de límites de variables
    v_min = [u_meas_min; x_min];
    v_max = [u_meas_max; x_max];
    
    % Límites de restricciones dinámicas (costura)
    tol = 1e-4;
    g_min = repmat(-tol, nx*N_e, 1);
    g_max = repmat(tol, nx*N_e, 1);
    
    % =====================================================================
    % CONFIGURACIÓN Y CONSTRUCCIÓN DEL SOLVER
    % =====================================================================
    
    % Vector de variables de decisión
    decvar = vertcat(u_est, x_est);
    
    % Vector de parámetros
    param = vertcat(y_meas_buff, u_buff, u_resist, x_est_last);
    
    % Definición del problema de optimización
    prob = struct('f', J, ...       % Función objetivo
                  'x', decvar, ...  % Variables de decisión
                  'g', g, ...       % Restricciones
                  'p', param);      % Parámetros
    
    % Configuración del solver IPOPT
    opts_solver.print_time = false;
    opts_solver.ipopt = struct(...
        'print_level', 0, ...                               % Sin salida por consola
        'linear_solver', 'ma27', ...                        % Solver lineal
        'fixed_variable_treatment', 'make_constraint', ...  % Tratamiento variables fijas
        'warm_start_init_point', 'yes', ...                 % Inicio caliente
        'tol', 1e-4, ...                                    % Tolerancia principal
        'acceptable_tol', 1e-4, ...                         % Tolerancia aceptable
        'max_wall_time', 2.5);                              % Tiempo máximo [s]
    
    % Crear solver
    solver = nlpsol('MHE_solver', 'ipopt', prob, opts_solver);
end
