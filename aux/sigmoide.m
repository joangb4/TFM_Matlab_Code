function y = sigmoide(x, x_c, sigma)
    % SIGMOIDE - Función sigmoidal desplazada y escalada
    %
    % Descripción:
    %   Calcula una función sigmoide (curva en forma de S) que transiciona
    %   suavemente entre 0 y 1
    %
    % Entradas:
    %   x     - Variable independiente (punto de evaluación)
    %   x_c   - Punto central de la sigmoide (donde y = 0.5)
    %   sigma - Factor de pendiente (controla la suavidad de la transición)
    %           Valores grandes de sigma -> transición más abrupta
    %           Valores pequeños de sigma -> transición más suave
    %
    % Salidas:
    %   y     - Valor de la sigmoide en x, en el intervalo (0, 1)
    %
    % Comportamiento:
    %   - Cuando x << x_c: y ~ 0
    %   - Cuando x = x_c:  y = 0.5
    %   - Cuando x >> x_c: y ~ 1

    y = 1/(1 + exp(-sigma*(x - x_c)));
end
