function [u_ctrl_seq, x_pred, lam_g, lam_x] = inicializar_controlador(N_c, size_g, x_k, u_k)
    % INICIALIZAR_CONTROLADOR - Inicialización del controlador MPC
    %
    % Descripción:
    %   Inicializa las variables del controlador con valores por defecto
    %   basados en el estado inicial y la primera acción de control
    %
    % Entradas:
    %   N_c        - Horizonte de control [p]
    %   size_g     - Tamaño del vector de restricciones
    %   x_k        - Estado inicial [nx x 1]
    %   u_k        - Primera entrada [nu x 1]
    %
    % Salidas:
    %   u_ctrl_seq - Secuencia de acciones de control iniciales
    %   x_pred     - Predicción de estado inicial para acción integral
    %   lam_g      - Multiplicadores de Lagrange (restricciones)
    %   lam_x      - Multiplicadores de Lagrange (límites)
    
    % Vector de secuencia de control a optimizar (replicar entrada inicial)
    u_ctrl_seq = repmat(u_k, N_c, 1);

    % Inicialización de estado predicho para acción integral
    x_pred = x_k;
    
    % Inicialización de multiplicadores de Lagrange
    lam_g = zeros(size_g, 1);
    lam_x = zeros(size(u_ctrl_seq, 1), 1);
end
