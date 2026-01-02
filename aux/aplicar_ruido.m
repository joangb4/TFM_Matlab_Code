function y_noisy = aplicar_ruido(y, sigmas)
    % APLICAR_RUIDO - Añade ruido gaussiano a una señal
    %
    % Descripción:
    %   Añade ruido aleatorio gaussiano (distribución normal) a una señal
    %   o conjunto de datos. Útil para simular mediciones con incertidumbre,
    %   ruido experimental, o perturbaciones en sistemas físicos
    %
    % Entradas:
    %   y      - Señal o datos originales (escalar, vector o matriz)
    %   sigmas - Desviación estándar del ruido gaussiano.
    %            Valores grandes de sigmas -> mayor amplitud de ruido
    %            Valores pequeños de sigmas -> ruido más sutil
    %
    % Salidas:
    %   y_noisy - Señal con ruido añadido, mismo tamaño que y
    %
    % Comportamiento:
    %   El ruido añadido sigue una distribución normal N(0, sigmas²)
    %   con media cero y desviación estándar sigmas. Cada elemento de
    %   la salida es: y_noisy(i) = y(i) + ruido_gaussiano(i)

    y_noisy = y + sigmas.*randn(size(y));
end
