function ode = modelo_mhe
    % MODELO_MHE - Modelo dinámico del sistema para estimador MHE
    %
    % Descripción:
    %   Implementa el modelo matemático del sistema. El modelo incluye:
    %   - Balance volumétrico (niveles de agua)
    %   - Balance energético (temperaturas)
    %
    % Estados (nx = 6):
    %   x(1) - h_w_C1:   Nivel agua tanque C1 [m]
    %   x(2) - UA:       Coef. transmisión calor [kW/°C]
    %   x(3) - d_UA:     Derivada de UA [kW/(°C·s)]
    %   x(4) - T_in_D1:  Temperatura entrada D1 [°C]
    %   x(5) - T_out_D1: Temperatura salida D1 [°C]
    %   x(6) - T_out_C1: Temperatura salida C1 [°C]
    %
    % Entradas (nu = 3):
    %   u(1) - Q_1_Lmin: Caudal 1 (bomba G1) [L/min]
    %   u(2) - a_v_2:    Apertura válvula 2 [%]
    %   u(3) - w_vent:   Velocidad ventilador [%]
    
    import casadi.*
    
    % =====================================================================
    % DEFINICIÓN DE PARÁMETROS Y VARIABLES SIMBÓLICAS
    % =====================================================================
    
    nx = 6; % Número de estados
    nu = 3; % Número de entradas
    
    % Variables simbólicas
    param_sim = SX.sym('param_sim', 10);    % Parámetros del modelo
    x         = SX.sym('x', nx);            % Vector de estados
    u         = SX.sym('u', nu);            % Vector de entradas
    u_resist  = SX.sym('u_resist', 1);      % Control resistencia [ON/OFF]
    
    % Asignación de parámetros del modelo
    r_D1     = param_sim(1);    % Radio tanque D1 [m]
    h_D1     = param_sim(2);    % Altura tanque D1 [m]
    r_C1     = param_sim(3);    % Radio tanque C1 [m]
    h_C1     = param_sim(4);    % Altura tanque C1 [m]
    rho_w    = param_sim(5);    % Densidad agua [kg/m³]
    Cp_w     = param_sim(6);    % Calor específico agua [kJ/(°C·kg)]
    T_amb    = param_sim(7);    % Temperatura ambiente [°C]
    M_w_refr = param_sim(8);    % Masa agua radiador [kg]
    M_w_tot  = param_sim(9);    % Masa agua total D1 + C1 [kg]
    W_resist = param_sim(10);   % Potencia resistencia [kW]
    
    % Asignación de estados
    h_w_C1   = x(1);    % Nivel agua C1 [m]
    UA       = x(2);    % Coef. transmisión calor [kW/°C]
    d_UA     = x(3);    % Derivada UA [kW/(°C·s)]
    T_in_D1  = x(4);    % Temperatura entrada D1 [°C]
    T_out_D1 = x(5);    % Temperatura salida D1 [°C]
    T_out_C1 = x(6);    % Temperatura salida C1 [°C]
    
    % Asignación de entradas de control
    Q_1_Lmin = u(1);      % Caudal bomba 1 [L/min]
    a_v_2    = u(2);      % Apertura válvula 2 [%]
    w_vent   = u(3) - 17; % Velocidad ventilador [%] (offset)
    
    % =====================================================================
    % MODELO VOLUMÉTRICO
    % =====================================================================

    % Área transversal de tanque C1
    A_C1 = pi*r_C1^2;   % [m²]
    A_D1 = pi*r_D1^2;   % [m²]

    % Masas de agua por tanque
    M_w_C1 = rho_w*h_w_C1*A_C1; % [kg]
    M_w_D1 = M_w_tot - M_w_C1;  % [kg]

    % Nivel de tanque D1
    h_w_D1 = M_w_D1/rho_w/A_D1; % [m]

    % Cálculo de caudal volumétrico 2 basado en apertura de válvula 2
    Q_2_Lmin = 0.0013*a_v_2^2 + 0.03*a_v_2; % [L/min]

    % Bloqueo de caudales cuando h_w = 0 (mediante sigmoide)
    h_min = 0.01;   % [m]
    slope = log(1/1e-5 - 1)/h_min;
    Q_1_Lmin = sigmoide(h_w_D1, h_min, slope)*Q_1_Lmin;
    Q_2_Lmin = sigmoide(h_w_C1, h_min, slope)*Q_2_Lmin;

    % Conversión de unidades a SI
    Q_1 = Q_1_Lmin/60/1000; % [m³/s]
    Q_2 = Q_2_Lmin/60/1000; % [m³/s]
    
    % Caudales másicos
    G_1 = rho_w*Q_1;    % [kg/s]
    G_2 = rho_w*Q_2;    % [kg/s]
    
    % Balance volumétrico en tanque C1
    d_h_w_C1 = (Q_1 - Q_2)/A_C1;    % [m/s]
    
    % =====================================================================
    % MODELO ENERGÉTICO - BALANCE DE ENERGÍA TÉRMICA
    % =====================================================================
    
    % Constante de tiempo para filtro dinámico de UA
    tau_UA = 2; % [s]
    
    % Modelo polinomial del coeficiente UA (obtenido experimentalmente)
    % UA = f(G_2, w_vent) - Correlación empírica
    UAs = 0.359*G_2 + 0.001919*w_vent - 5.285*G_2^2 + ...
          0.16*G_2*w_vent - 0.0001674*w_vent^2 + 30.31*G_2^3 - ...
          0.882*G_2^2*w_vent - 0.002611*G_2*w_vent^2 + ...
          5.029e-06*w_vent^3 - 56.23*G_2^4 + 1.761*G_2^3*w_vent + ...
          0.009239*G_2^2*w_vent^2 + 2.037e-05*G_2*w_vent^3 - ...
          6.364e-08*w_vent^4 - 0.6429*G_2^4*w_vent - ...
          0.01134*G_2^3*w_vent^2 - 2.92e-05*G_2^2*w_vent^3 - ...
          6.462e-08*G_2*w_vent^4 + 2.888e-10*w_vent^5;
    
    % Implementación de filtro de segundo orden para UA (añadir dinámica)
    d2_UA = (UAs - UA)/tau_UA - 2*tau_UA/tau_UA*d_UA;   % [kW/(°C·s)]
    
    % Cálculo de diferencias de temperatura para DTML
    diff_T_out_C1 = T_out_C1 - T_amb;   % [°C]
    diff_T_in_D1  = T_in_D1 - T_amb;    % [°C]

    % Regularización para evitar problemas numéricos (softplus)
    epsilon_T = 1e-6;
    diff_T_out_C1 = sqrt(diff_T_out_C1^2 + epsilon_T);
    diff_T_in_D1  = sqrt(diff_T_in_D1^2 + epsilon_T);
        
    % Diferencia de temperatura media logarítmica (DTML) - Aproximación Chen
    DTML = (diff_T_out_C1*diff_T_in_D1*(diff_T_out_C1 + diff_T_in_D1)/2)^(1/3); % [°C]
    
    % Potencia de refrigeración del radiador
    W_refr = UA*DTML;   % [kW]
    
    % Balance energético
    % Radiador (intercambiador de calor)
    d_T_in_D1 = ((T_out_C1 - T_in_D1)*G_2*Cp_w - W_refr)/(Cp_w*M_w_refr);   % [°C/s]
    
    % Tanque D1 (con resistencia calefactora)
    d_T_out_D1 = (G_2*(T_in_D1 - T_out_D1)*Cp_w + u_resist*W_resist)/(M_w_D1*Cp_w); % [°C/s]
    
    % Tanque C1
    d_T_out_C1 = G_1*(T_out_D1 - T_out_C1)/M_w_C1;  % [°C/s]
    
    % =====================================================================
    % ENSAMBLE DEL SISTEMA DE ECUACIONES DIFERENCIALES
    % =====================================================================
    
    dxdt    = SX.zeros(nx, 1);
    dxdt(1) = d_h_w_C1;    % Derivada nivel C1
    dxdt(2) = d_UA;        % Derivada UA
    dxdt(3) = d2_UA;       % Derivada segunda UA
    dxdt(4) = d_T_in_D1;   % Derivada temperatura entrada D1
    dxdt(5) = d_T_out_D1;  % Derivada temperatura salida D1
    dxdt(6) = d_T_out_C1;  % Derivada temperatura salida C1

    % Parámetros de estructura
    param = vertcat(u, u_resist, param_sim);
    
    % Estructura de salida para CasADi
    ode = struct('x', x, ...
                 'p', param, ...
                 'ode', dxdt);
end
