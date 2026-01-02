function [v_est_seq, x_est, u_est, lam_g_est, lam_x_est, y_meas_buff, u_buff, status, n_iter] = estimar_estado(N_e, solver, v_est_seq_0, lam_g_est_0, lam_x_est_0, v_est_min, v_est_max, g_est_min, g_est_max, y_meas_buff, y_meas_k, u_buff, u_k, u_resist)
    import casadi.*

    % ESTIMAR_ESTADO
    %
    % Descripción:
    %   Ejecuta la estimación para un instante de tiempo, actualizando
    %   los históricos y resolviendo el problema de estimación.
    %
    % Entradas:
    %   N_e           - Horizonte de estimación
    %   solver        - Solver IPOPT configurado
    %   v_est_seq_0   - Secuencia de variables estimadas inicial [(nx + 2)*N_e x 1]
    %   lam_g_est_0   - Multiplicadores de restricciones
    %   lam_x_est_0   - Multiplicadores de límites
    %   v_est_min/max - Límites de variables
    %   g_est_min/max - Límites de restricciones
    %   y_meas_buff   - Histórico de salidas medidas [3*N_e x 1]
    %   y_meas_k      - Nueva medición [3 x 1]
    %   u_buff        - Histórico de entradas [nu*N_e x 1]
    %   u_k           - Nueva entrada [nu x 1]
    %   u_resist      - Control resistencia [ON/OFF]
    %
    % Salidas:
    %   v_est_seq     - Secuencia de variables estimadas óptima
    %   x_est         - Estados estimados en instante actual
    %   u_est         - Entradas estimadas [Q_1, a_v_2]
    %   lam_g_est     - Multiplicadores actualizados de restricciones
    %   lam_x_est     - Multiplicadores actualizados de límites
    %   y_meas_buff   - Histórico de salidas medidas actualizado
    %   u_buff        - Histórico de entradas actualizado
    %   status        - Estado de convergencia del solver
    %   n_iter        - Número de iteraciones de estimación
    
    % =====================================================================
    % DEFINICIÓN DE PARÁMETROS DEL SISTEMA
    % =====================================================================
    
    nx = 6; % Número de estados
    nu = 3; % Número de entradas
    
    % =====================================================================
    % ACTUALIZACIÓN DE HISTÓRICOS CON NUEVAS MEDICIONES
    % =====================================================================
    
    % Actualizar históricos de salidas medidas (horizonte móvil)
    y_meas_buff = [y_meas_buff(4:end); y_meas_k];
    
    % Actualizar históricos de entradas (horizonte móvil)
    u_buff = [u_buff(nu + 1 : end); u_k];
    
    % =====================================================================
    % CONSTRUCCIÓN DEL VECTOR DE PARÁMETROS
    % =====================================================================
    
    % Extraer estimación anterior de estados (para regularización)
    x_est_0 = v_est_seq_0(2*N_e + 1 : end);
    
    % Ensamblar vector de parámetros para el solver
    param_est = vertcat(y_meas_buff, ...    % Histórico de salidas medidas
                        u_buff, ...         % Histórico de entradas
                        u_resist, ...       % Control resistencia
                        x_est_0);           % Estimación anterior
    
    % =====================================================================
    % ENSAMBLE DEL SOLVER DE OPTIMIZACIÓN
    % =====================================================================
    
    sol = solver('x0', v_est_seq_0, ...     % Secuencia de variables estimadas inicial
                 'lam_g0', lam_g_est_0, ... % Multiplicadores restricciones
                 'lam_x0', lam_x_est_0, ... % Multiplicadores límites
                 'p', param_est, ...        % Parámetros del solver
                 'lbx', v_est_min, ...      % Límites inferiores variables
                 'ubx', v_est_max, ...      % Límites superiores variables
                 'lbg', g_est_min, ...      % Límites inferiores restricciones
                 'ubg', g_est_max);         % Límites superiores restricciones
        
    % =====================================================================
    % EXTRACCIÓN DE RESULTADOS DEL SOLVER
    % =====================================================================
    
    % Extraer solución óptima
    v_est_seq = full(sol.x);
    lam_g_est = full(sol.lam_g);
    lam_x_est = full(sol.lam_x);

    % Extraer información del solver
    status         = solver.stats.return_status;
    n_iter         = solver.stats.iter_count;
    solver_success = solver.stats.success;

    % Extraer soluciones de la secuencia
    x_est = v_est_seq(end - nx + 1 : end);              % Estado estimado
    u_est = [v_est_seq(2*N_e - 1); v_est_seq(2*N_e)];   % Entradas estimadas

    % =====================================================================
    % ESTRATEGIA DE RECUPERACIÓN EN CASO DE FALLO
    % =====================================================================
    
    if ~ solver_success 
        fprintf('MHE: Fallo en optimización principal. Estado: %s\n', status);
        
        x_est(1) = y_meas_k(1);
        x_est(4) = y_meas_k(2);
        x_est(5) = y_meas_k(3);
    end    
end
