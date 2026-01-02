function integrador = constructor_integrador(ode, dt)
    % CONSTRUCTOR_INTEGRADOR - Crea integrador numérico para simulación
    %
    % Descripción:
    %   Construye un integrador CVODES para la simulación del modelo dinámico
    %
    % Entradas:
    %   ode        - Estructura con el modelo dinámico (ecuaciones diferenciales)
    %   dt         - Paso de integración [s]
    %
    % Salidas:
    %   integrador - Objeto integrador de CasADi configurado
    
    import casadi.*
    
    % Configuración de tolerancias para el integrador CVODES
    opts = struct('abstol', 1e-4, ...  % Tolerancia absoluta
                  'reltol', 1e-4);     % Tolerancia relativa
    
    % Crear integrador con método CVODES
    integrador = integrator('F', 'cvodes', ode, 0, dt, opts);
end
