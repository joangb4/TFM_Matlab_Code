function x_kp1 = simulador_1p(simulador, parametros, x_k, u_k, u_resist)
    % SIMULADOR_1P - Ejecuta un paso de simulación del modelo dinámico
    %
    % Descripción:
    %   Avanza el estado del sistema un intervalo temporal dt utilizando
    %   el integrador numérico CVODES previamente configurado
    %
    % Entradas:
    %   simulador  - Objeto integrador de CasADi
    %   parametros - Vector de parámetros del modelo
    %   x_k        - Estado actual del sistema en el instante k
    %   u_k        - Entrada de control en el instante k
    %   u_resist   - Entrada de resistencia
    %
    % Salidas:
    %   x_kp1      - Estado del sistema en el instante k + 1 (tras avanzar dt)
    
    % Parámetros de simulación
    param_sim = vertcat(u_k, u_resist, parametros);

    % Ejecutar integración numérica desde t = 0 hasta t = dt
    sol = simulador('x0', x_k, 'p', param_sim);

    % Estado final
    x_kp1 = full(sol.xf);
end
