function [solver, u_min, u_max, g_min, g_max] = constructor_mpc_SS(N_h, N_c, W_delta, param_modelo, simulador)
    % CONSTRUCTOR_MPC_SS - Construye optimizador MPC con Single Shooting
    %
    % Descripción:
    %   Implementa la formulación de Model Predictive Control (MPC) usando
    %   la técnica de Single Shooting para optimización de acciones de
    %   control óptimas.
    %
    % Entradas:
    %   N_h          - Horizonte de predicción [pa]
    %   N_c          - Horizonte de control [p]
    %   R_ctrl       - Factor de penalización de control frente a seguimiento
    %   W_track      - Matriz de pesos de seguimiento de referencias
    %   W_delta      - Matriz de pesos de penalización de cambios en el control
    %   param_modelo - Parámetros del modelo físico
    %   simulador    - Integrador numérico del modelo
    %
    % Salidas:
    %   solver       - Optimizador IPOPT configurado
    %   u_min        - Límites inferiores de acciones de control
    %   u_max        - Límites superiores de acciones de control
    %   g_min        - Límites inferiores de restricciones
    %   g_max        - Límites superiores de restricciones

    import casadi.*

    % =====================================================================
    % DEFINICIÓN DE PARÁMETROS DEL PROBLEMA
    % =====================================================================
    
    nx = 6; % Número de estados
    nu = 3; % Número de entradas

    % Extracción de parámetros físicos del modelo
    p_r_D1  = param_modelo(1);  % Radio tanque D1 [m]
    p_h_D1  = param_modelo(2);  % Altura tanque D1 [m]
    p_r_C1  = param_modelo(3);  % Radio tanque C1 [m]
    p_h_C1  = param_modelo(4);  % Altura tanque C1 [m]
    p_rho_w = param_modelo(5);  % Densidad agua [kg/m³]
    p_T_amb = param_modelo(7);  % Temperatura ambiente [°C]
    p_A_C1  = pi*p_r_C1^2;      % Área transversal de tanque C1 [m²]
    p_A_D1  = pi*p_r_D1^2;      % Área transversal de tanque D1 [m²]

    % =====================================================================
    % VARIABLES DE DECISIÓN DEL PROBLEMA MPC
    % =====================================================================
    
    % Variables a estimar
    u_ctrl = MX.sym('u_ctrl', nu*N_c);  % Acciones de control [a_v_1, a_v_2, w_vent]
    
    % Parámetros del problema (datos conocidos)
    x_ini    = MX.sym('x_ini', nx);     % Estado inicial
    u_ini    = MX.sym('u_ini', nu);     % Último control aplicado
    u_resist = MX.sym('u_resist', 1);   % Control resistencia [ON/OFF]
    ref      = MX.sym('ref', 2*N_h);    % Referencias [h_D1, T_out]
    x_pred   = MX.sym('x_pred', nx);    % Predicción anterior del modelo MPC
    w_track  = MX.sym('w_track', 2);    % Vector de pesos por error de seguimiento
    R_ctrl   = MX.sym('R_ctrl', 1);     % Peso global vs W_track
    
    % =====================================================================
    % INICIALIZACIÓN DE FUNCIÓN OBJETIVO Y RESTRICCIONES
    % =====================================================================
    
    J = 0;  % Función objetivo
    g = []; % Vector de restricciones

    % Matriz de pesos por error de seguimiento (h_w_D1 y T_out_D1) (normalizado)
    W_track = diag(w_track'./[p_h_D1, 80]);
    %              [h_w_D1 T_out_D1]
    
    % =====================================================================
    % ACCIÓN INTEGRAL - CORRECCIÓN DE ESTADO INICIAL
    % =====================================================================

    % Cálculo de error entre estimación MHE y predicción MPC anterior
    e_act = x_ini - x_pred;
    
    % Corrección de estado inicial con acción integral
    x_k = x_ini + e_act;

    u_km1 = u_ini;  % Cambio de nomenclatura a u_{k-1}

    % =====================================================================
    % BUCLE DE HORIZONTE DE CONTROL (N_c)
    % =====================================================================

    for k = 1:N_c
        
        % Indexación de variables en instante k
        idx_u_k   = (k - 1)*nu + 1 : k*nu;  % Acción de control u_k
        idx_ref_k = (k - 1)*2 + 1 : k*2;    % Referencias [h_w_D1, T_out_D1]
        
        % Extracción de variables
        u_k   = u_ctrl(idx_u_k);    % Control a optimizar
        ref_k = ref(idx_ref_k);     % Referencia actual
        
        % Propagación del modelo dinámico
        param_sim_k = vertcat(u_k, u_resist, param_modelo);
        sim_res     = simulador('x0', x_k, 'p', param_sim_k);
        x_kp1       = sim_res.xf;   % Estado predicho x_{k+1}

        % Conversión de h_w_C1 a h_w_D1
        M_w_C1 = p_rho_w*x_kp1(1)*p_A_C1;   % Masa de agua tanque C1 [kg]
        M_w_D1 = 13.7 - M_w_C1;             % Masa de agua tanque D1 [kg]
        h_w_D1 = M_w_D1/p_rho_w/p_A_D1;     % Nivel de tanque D1 [m]
        
        % Construcción de salidas controladas
        y_ctrl_k = [h_w_D1;     % h_w_D1
                    x_kp1(5)];  % T_out_D1
        
        % Término de seguimiento de referencias en función objetivo
        e = [y_ctrl_k(1) - ref_k(1);    % Error cuadrático en altura
             y_ctrl_k(2) - ref_k(2)];       % Error lineal en temperatura
        J = J + e'*W_track*e;
        
        % Término de penalización de variación de control en función objetivo
        du = u_km1 - u_k;                               % Cambio en acción de control
        J  = J + R_ctrl*(du'*W_delta*du + 0.05*u_k(3)); % Penalización adicional del ventilador
        
        % Restricciones de salida
        g = [g; h_w_D1;     % h_w_D1
                x_kp1(1);   % h_w_C1
                x_kp1(4);   % T_in_D1
                x_kp1(5);   % T_out_D1
                x_kp1(6)];  % T_out_C1
        
        % Actualización para siguiente iteración
        u_km1 = u_k;    % Control
        x_k   = x_kp1 + e_act;  % Estados
    end
    
    % =====================================================================
    % BUCLE DE HORIZONTE DE PREDICCIÓN (SIN CONTROL) (N_c + 1 : N_h)
    % =====================================================================
    
    for k = (N_c + 1):N_h
        
        % Indexación de referencias
        idx_ref_k = (k - 1)*2 + 1 : k*2;
        
        % Extracción de variables
        ref_k = ref(idx_ref_k);
        
        % Propagación del modelo dinámico con última acción de control
        param_sim_k = vertcat(u_k, u_resist, param_modelo);
        sim_res     = simulador('x0', x_k, 'p', param_sim_k);
        x_kp1       = sim_res.xf;
        
        % Conversión de h_w_C1 a h_w_D1
        M_w_C1 = p_rho_w*x_kp1(1)*p_A_C1;   % Masa de agua tanque C1 [kg]
        M_w_D1 = 13.7 - M_w_C1;             % Masa de agua tanque D1 [kg]
        h_w_D1 = M_w_D1/p_rho_w/p_A_D1;     % Nivel de tanque D1 [m]
        
        % Construcción de salidas controladas
        y_ctrl_k = [h_w_D1;     % h_w_D1
                    x_kp1(5)];  % T_out_D1
        
        % Término de seguimiento de referencias en función objetivo
        e = [y_ctrl_k(1) - ref_k(1);        % Error cuadrático en altura
             y_ctrl_k(2) - ref_k(2)];       % Error lineal en temperatura
        J = J + e'*W_track*e;
        
        % Restricciones de salida
        g = [g; h_w_D1;     % h_w_D1
                x_kp1(1);   % h_w_C1
                x_kp1(4);   % T_in_D1
                x_kp1(5);   % T_out_D1
                x_kp1(6)];  % T_out_C1
        
        % Actualización para siguiente iteración
        x_k = x_kp1 + e_act;    % Estados
    end
    
    % =====================================================================
    % DEFINICIÓN DE LÍMITES Y RESTRICCIONES
    % =====================================================================
    
    % Límites de acciones de control
    u_min = repmat([20;         % a_v_1 mínimo
                    20;         % a_v_2 mínimo
                    17], ...    % w_vent mínimo
                    N_c, 1);
    
    u_max = repmat([75;         % a_v_1 máximo
                    100;        % a_v_2 máximo
                    100], ...   % w_vent máximo
                    N_c, 1);
    
    % Límites de restricciones de salida
    g_min = repmat([0.1*p_h_D1;         % h_w_C1 mínimo
                    0.1*p_h_C1;         % h_w_C1 mínimo
                    p_T_amb - 5;        % T_in_D1 mínimo
                    p_T_amb - 5;        % T_out_D1 mínimo
                    p_T_amb - 5], ...   % T_out_C1 mínimo
                    N_h, 1);
    
    g_max = repmat([0.95*p_h_D1;    % h_w_C1 máximo
                    0.92*p_h_C1;    % h_w_C1 máximo
                    60;             % T_in_D1 máximo
                    60;             % T_out_D1 máximo
                    60], ...        % T_out_C1 máximo
                    N_h, 1);
    
    % =====================================================================
    % CONFIGURACIÓN Y CONSTRUCCIÓN DEL SOLVER
    % =====================================================================
    
    % Vector de variables de decisión
    decvar = u_ctrl;
    
    % Vector de parámetros
    param = vertcat(x_ini, u_ini, u_resist, ref, x_pred, w_track, R_ctrl);
    
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
    solver = nlpsol('MPC_solver', 'ipopt', prob, opts_solver);
end
