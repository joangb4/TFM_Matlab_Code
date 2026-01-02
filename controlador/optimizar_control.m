function [u_ctrl_seq, u_ctrl_k, x_pred, lam_g_ctrl, lam_x_ctrl, status, n_iter] = optimizar_control(solver, u_ctrl_seq_0, lam_g_ctrl_0, lam_x_ctrl_0, u_ctrl_min, u_ctrl_max, g_ctrl_min, g_ctrl_max, x_est_k, u_k, u_resist, ref_k, x_pred_0, y_meas_k, simulador, w_track, R_ctrl, param_modelo)
    import casadi.*
    
    % OPTIMIZAR_CONTROL
    %
    % Descripción:
    %   Ejecuta el solver IPOPT para obtener la secuencia óptima de control
    %   que minimiza la función objetivo del MPC sujeta a restricciones.
    %   Incluye estrategia de contingencia en caso de fallo del optimizador.
    %
    % Entradas:
    %   solver         - Solver IPOPT configurado
    %   u_ctrl_seq_0   - Secuencia de control inicial [nu*N_c x 1]
    %   lam_g_ctrl_0   - Multiplicadores de restricciones
    %   lam_x_ctrl_0   - Multiplicadores de límites
    %   u_ctrl_min/max - Límites de acción de control
    %   g_ctrl_min/max - Límites de restricciones
    %   x_est_k        - Estado estimado actual [nx x 1]
    %   u_k            - Última acción de control aplicada [nu x 1]
    %   u_resist       - Control resistencia [ON/OFF]
    %   ref_k          - Vector de referencias [2*N_h x 1]
    %   x_pred_0       - Predicción anterior del estado [nx x 1]
    %   y_meas_k       - Mediciones actuales [3 x 1]
    %   simulador      - Integrador numérico del modelo
    %   param_modelo   - Parámetros del modelo físico
    %
    % Salidas:
    %   u_ctrl_seq     - Secuencia de control óptima [nu*N_c x 1]
    %   u_ctrl_k       - Acción de control óptimo a aplicar [nu x 1]
    %   x_pred         - Predicción de estado para siguiente iteración [nx x 1]
    %   lam_g_ctrl     - Multiplicadores actualizados de restricciones
    %   lam_x_ctrl     - Multiplicadores actualizados de límites
    %   status         - Estado de convergencia del solver
    %   n_iter         - Número de iteraciones de optimización
    
    % =====================================================================
    % DEFINICIÓN DE PARÁMETROS DEL SISTEMA
    % =====================================================================
    
    nu = 3; % Número de entradas
    
    % =====================================================================
    % CONSTRUCCIÓN DEL VECTOR DE PARÁMETROS
    % =====================================================================
    
    % Ensamble de parámetros para el solver MPC
    param_ctrl = vertcat(x_est_k, ...   % Estado inicial estimado
                         u_k, ...       % Control previo
                         u_resist, ...  % Control resistencia
                         ref_k, ...     % Referencias
                         x_pred_0, ...  % Predicción anterior
                         w_track', ...  % Vector de penalización de seguimiento
                         R_ctrl);       % Peso global vs W_track
    % =====================================================================
    % ENSAMBLE DEL SOLVER DE OPTIMIZACIÓN
    % =====================================================================
    
    sol = solver('x0', u_ctrl_seq_0, ...        % Secuencia de control inicial
                 'lam_g0', lam_g_ctrl_0, ...    % Multiplicadores restricciones
                 'lam_x0', lam_x_ctrl_0, ...    % Multiplicadores límites
                 'p', param_ctrl, ...           % Parámetros del solver
                 'lbx', u_ctrl_min, ...         % Límites inferiores control
                 'ubx', u_ctrl_max, ...         % Límites superiores control
                 'lbg', g_ctrl_min, ...         % Límites inferiores restricciones
                 'ubg', g_ctrl_max);            % Límites superiores restricciones
    
    % =====================================================================
    % EXTRACCIÓN DE RESULTADOS DEL SOLVER
    % =====================================================================
    
    % Extraer multiplicadores óptimos
    lam_g_ctrl = full(sol.lam_g);
    lam_x_ctrl = full(sol.lam_x);
    
    % Extraer información del solver
    status         = solver.stats.return_status;
    n_iter         = solver.stats.iter_count;
    solver_success = solver.stats.success;
    
    % =====================================================================
    % ESTRATEGIA DE RECUPERACIÓN EN CASO DE FALLO
    % =====================================================================
    
    if ~ solver_success
        % El solver no convergió: aplicar estrategia de emergencia
        u_ctrl_seq = u_ctrl_seq_0;
        
        % Estrategia basada en mediciones
        if y_meas_k(1) > g_ctrl_max(1)
            % Reducir entrada y aumentar salida de agua en tanque C1
            u_ctrl_k = [20;                     % a_v_1 mínimo
                        75;                     % a_v_2 máximo
                        u_ctrl_seq(3 + nu)];    % Mantener velocidad ventilador
            
        elseif y_meas_k(1) < g_ctrl_min(1)
            % Aumentar entrada y reducir salida de agua en tanque C1
            u_ctrl_k = [75;                     % a_v_1 máximo
                        20;                     % a_v_2 mínimo
                        u_ctrl_seq(3 + nu)];    % Mantener velocidad ventilador
        else
            % Mantener control anterior
            u_ctrl_k = [u_ctrl_seq(1 + nu);     % a_v_1 predicho
                        u_ctrl_seq(2 + nu);     % a_v_2 predicho
                        u_ctrl_seq(3 + nu)];    % w_vent predicho
        end
        
    else
        % =================================================================
        % EXTRACCIÓN DE SOLUCIÓN ÓPTIMA
        % =================================================================
        
        % Extraer solución óptima
        u_ctrl_seq = full(sol.x);
        
        % Extraer solución de la secuencia
        u_ctrl_k = [u_ctrl_seq(1);  % a_v_1
                    u_ctrl_seq(2);  % a_v_2
                    u_ctrl_seq(3)]; % vel_vent
    end
    
    % =====================================================================
    % PREDICCIÓN DE ESTADO PARA ACCIÓN INTEGRAL
    % =====================================================================
    
    % Parámetros de simulación con control a aplicar
    param_sim = vertcat(u_ctrl_k, u_resist, param_modelo);
    
    % Simulación de un paso con el control óptimo
    sol_pred = simulador('x0', x_est_k, 'p', param_sim);
    
    % Estado predicho para siguiente iteración (acción integral MPC - MHE)
    x_pred = full(sol_pred.xf);
end
