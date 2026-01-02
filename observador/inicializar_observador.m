function [v_est_seq, y_meas_buff, u_buff, lam_g, lam_x] = inicializar_observador(N_e, size_g, x_k, y_meas_k, u_k)
    % INICIALIZAR_OBSERVADOR - Inicialización del estimador MHE
    %
    % Descripción:
    %   Inicializa las variables del observador con valores por defecto
    %   basados en el estado inicial y las primeras mediciones
    %
    % Entradas:
    %   N_e         - Horizonte de estimación [p]
    %   size_g      - Tamaño del vector de restricciones
    %   x_k         - Estado inicial [nx x 1]
    %   y_meas_k    - Primera medición [3 x 1]
    %   u_k         - Primera entrada [nu x 1]
    %
    % Salidas:
    %   v_est_seq   - Secuencia de variables estimadas iniciales
    %   y_meas_buff - Histórico de salidas medidas inicial
    %   u_buff      - Histórico de entradas inicial
    %   lam_g       - Multiplicadores de Lagrange (restricciones)
    %   lam_x       - Multiplicadores de Lagrange (límites)
    
    % Inicialización de estados (replicar estado inicial)
    x_est = repmat(x_k, (N_e + 1), 1);
    
    % Inicialización de entradas estimadas [Q_1_Lmin, a_v_2]
    u_est_init = repmat(u_k(1:2), N_e, 1);
    
    % Vector de secuencia de variables a estimar
    v_est_seq = vertcat(u_est_init, x_est);
    
    % Inicialización de buffers
    y_meas_buff = repmat(y_meas_k, N_e, 1);
    u_buff      = repmat(u_k, N_e, 1);
    
    % Inicialización de multiplicadores de Lagrange
    lam_g = zeros(size_g, 1);
    lam_x = zeros(size(v_est_seq, 1), 1);
end
