%Columbia University Rocketry
%Ideal Nitrous Model (Margeret Mary's Masters Thesis
%Burnell's choked flow model for nitrous injection
%Explicit Method - Forward Difference

% TODO: Change units to imperial
% TODO: Replace interpolation of C* and Thrust Coeff with equations

clear all
clc

global OUTPUT_RUN_GRAPHS;

RUN_ONCE = 1;
OUTPUT_RUN_GRAPHS = 1;
OUTPUT_GRADIENT_GRAPHS = 1;
OUTPUT_IMPULSE_GRAPH = 1;

OUTPUT_THRUST_CURVE = 1;

%******************************THRUST CURVE DATA************************************%
% These are not used in these calculations but are used in
% OpenRocket/RASAero simulations
MOTOR_NAME = "IREC_Bald_Guy_From_Die_Hard";
MOTOR_DIAMETER = "146.05";      %Diameter of motor (in mm)
MOTOR_LENGTH = "1651.0";        %Length of motor (in mm)
MOTOR_DELAYS = "27";            %Delay between ignition and parachute
MOTOR_MASS = 8.0;               %Mass of motor (kg)
MOTOR_PROP_WEIGHT = "7.82";     %Total Propellant weight (kg)
MOTOR_TOTAL_WEIGHT = "7.83";    %Total weight of motor+propellant (kg)
MOTOR_MANUFACTURER = "Elongated_Muskrat";

% Where to save thrust curve
MOTOR_FILEPATH = "C:\Users\dkola\AppData\Roaming\OpenRocket\ThrustCurves";
% What to name the file
MOTOR_FILENAME = "Bald_Guy_From_Die_Hard.eng";
%****************************END THRUST CURVE DATA**********************************%

MOTOR_FP = strcat(MOTOR_FILEPATH,"\",MOTOR_FILENAME);   %Create thrust curve file pointer

global DEnd;
global m_T;
global Fuel_Density;

DEnd = 0.08255;                 % Final Port Diameter [m]
m_T = 2.2882;                   % tank mass [kg]
Fuel_Density = 812.4;           % [kg/m^3]

global Thrust;

% Given constants

%****************************************************************
global m_loaded;
m_loaded = 6.7;    % N2O mass initially loaded into tank [kg]
Main_m_loaded = m_loaded;
%****************************************************************
% Might want to calculate total area based on injector hole diameter and
% number of holes just to make it easier to test different values
global Ainj;
Ainj = 0.0000535;    % injector area (sum of individual hole areas) [m^2]
Main_Ainj = Ainj;
%****************************************************************
global V;
V = 0.01;             % total tank volume [m^3]
Main_V = V;
%****************************************************************
Nozzle_Throat_Diameter = 0.025;                %Diameter of nozzle throat [m]
global A_Star;
A_Star = pi*(Nozzle_Throat_Diameter/2)*(Nozzle_Throat_Diameter/2);      % Nozzle Throat Area [m^3]
Main_A_Star = A_Star;
%****************************************************************
global Length;
Length = 0.24;           % Grain Length [m]
Main_Length = Length;
%****************************************************************
global D0;
D0 = 0.03175;              % Initial Port Diameter [m]
D = D0;

%Read data from RPA csv
%Make sure the file is in the same directory as you are
rawdata = csvread('nestedanal5.csv',1,0);

%The values used in the data from RPA, low, high and step values for OF and
%pressure
OFLow = 1.0;
OFHigh = 35.0;
OFStep = 0.1;

PLow = 101325;
PHigh = 9500000;
PStep = 100000;

global data_temp;
global data_gamma;
global data_Cstar;
global data_thrust;
global data_gridx;
global data_gridy;

      %[    OF        Chamber P    Chamber T      Gamma         C*      Cf]
data = [rawdata(:,1) rawdata(:,2) rawdata(:,6) rawdata(:,8) rawdata(:,10) rawdata(:,13)];
y = PLow:PStep:PHigh;
y = [y PHigh];
% Create the data structures used for the interpolation
% The ones for C* and Thrust Coeff are not really needed since these can be
% calculated using other known parameters. This will use much less
% computing power and be more accurate
% TODO: Replace interpolation of C* and Thrust Coeff with equations
[data_gridx,data_gridy] = meshgrid(OFLow:OFStep:OFHigh,y);
data_temp = reshape(data(:,3),[],1 + (OFHigh - OFLow)/OFStep);
data_gamma = reshape(data(:,4),[],1 + (OFHigh - OFLow)/OFStep);
data_Cstar = reshape(data(:,5),[],1 + (OFHigh - OFLow)/OFStep);     %Shouldn't need this
data_thrust = reshape(data(:,6),[],1 + (OFHigh - OFLow)/OFStep);    %Shouldn't need this

% Run the simulation and get the total impulse
totalImpulse = calcImpulse()

%Output thrust curve in .eng format
if OUTPUT_THRUST_CURVE == 1
    %calculate total weight of propellant
    %MOTOR_PROP_WEIGHT = num2str(m_loaded + Fuel_Density*pi*realpow((DEnd - D0)/2,2)*Length)
    %MOTOR_TOTAL_WEIGHT = num2str(MOTOR_MASS + m_loaded + Fuel_Density*pi*realpow((DEnd - D0)/2,2)*Length)
    %MOTOR_PROP_WEIGHT = num2str(7.82)
    %MOTOR_TOTAL_WEIGHT = num2str(0.01)
    %Create array of header values
    TC_Arr = [MOTOR_NAME MOTOR_DIAMETER MOTOR_LENGTH MOTOR_DELAYS MOTOR_PROP_WEIGHT ...
        MOTOR_TOTAL_WEIGHT MOTOR_MANUFACTURER];
    %Create array of thrust values
    T = [Thrust(:,1)'; Thrust(:,2)'];
    %Remove first row since this is a 0 row and that signals the end of a
    %.eng file
    T = T(:,2:end);
    TCFile = fopen(MOTOR_FP,'w');
    %output data to .eng file
    fprintf(TCFile,'; Thrust Curve Output by Chuckleton\n');
    fprintf(TCFile,'%s %s %s %s %s %s %s\n',TC_Arr);
    fprintf(TCFile,'   %f %f\n',T);
    fclose(TCFile);
end

%function to calculate the impulse given by the engine given a small change
%to one or more parameters from its original value
function imp = calcImpulse()
    Ti = 294.7;             % Initial nitrous bulk temperature
    R = 8314.3;             % universal gas constant [J/(kmol*K)]

    a_coeff = 0.00016;     % aG^n regression rate coeffs
    n_coeff = 0.5;         % Determined experimentally

    % Efficiency Coefficients
    Nozzle_Efficiency = 0.9772;
    Reaction_Efficiency = 0.77;
    Drag_Efficiency = 0.96223;

    %Injector Discharge Coefficient
    Cd = 0.425;

    P0 = 101325;             % Initial Pressure [Pa]
    Pe = P0;                % Chamber Pressure [Pa]
    
    MW2 = 44.013;           % molecular weight of N2O
    M = 26.2448;            % M of mixture in chamber [g/mol]
    R1 = R / M;             % Adjusted R coefficient
    
    n_He = 0;               % helium gas [kmol]

    
    global m_T;
    global Fuel_Density;
    global m_loaded;
    global Ainj;
    global V;
    global A_Star;
    global Length;
    global D0;
    D = D0;

    %Data used for interpolation
    global data_temp;
    global data_gamma;
    global data_Cstar;
    global data_thrust;
    global data_gridx;
    global data_gridy;
    
    global Thrust;
    
    % Perry's Chemical Engineers' Handbook Property Equations
    G1 = 96.512;        % vapour pressure of N2O [Pa] coefficients
    G2 = -4045;         % valid for Temp range [182.3 K - 309.57 K]
    G3 = -12.277;      
    G4 = 2.886e-5;
    G5 = 2;

    Tc = 309.57;        % critical temperature of N2O [K]
    J1 = 2.3215e7;      % heat of vapourization of N2O [J/kmol] coefficients
    J2 = 0.384;         % valid for Temp range [182.3K - 309.57 K]
    J3 = 0;
    J4 = 0;

    C1 = 0.2079e5;      % heat capacity of He at constant pressure [K/(kmol*K)] coefficients
    C2 = 0;             % valid for Temp range [100 K - 1500 K]
    C3 = 0;
    C4 = 0;
    C5 = 0;

    D1 = 0.2934e5;      % heat capacity of N2O gas at constant pressure [J/(kmol*K)] coefficients
    D2 = 0.3236e5;      % valid for Temp range [100 K - 200 K]
    D3 = 1.1238e3;      
    D4 = 0.2177e5;
    D5 = 479.4;

    E1 = 6.7556e4;      % heat capacity of N2O liquid at constant pressure [J/(kmol*K)] coefficients
    E2 = 5.4373e1;      % valid for Temp range [182.3 K - 200K]
    E3 = 0;
    E4 = 0;
    E5 = 0;

    Q1 = 2.781;         %molar specific volume of liquid N2O [m^3/kmol] coefficients
    Q2 = 0.27244;
    Q3 = 309.57;
    Q4 = 0.2882;

    % Initial conditions
    n_to = m_loaded/MW2;                                        % initial total N2O in tank [kmol]
    Vhat_li = Q2^(1+(1-Ti/Q3)^Q4)/Q1;                           % molar volume of liquid N2O [m^3/kmol]    
    To = Ti;                                                    % initial temperature [K]   
    P_sato = exp(G1 + G2/To + G3*log(To) + G4*To^G5);           % initial vapour pressure of N2O [Pa]
    n_go = P_sato*(V - Vhat_li*n_to)/(-P_sato*Vhat_li + R*To);  % initial N2O gas [kmol]
    n_lo = (n_to*R*To - P_sato*V)/(-P_sato*Vhat_li + R*To);     % initial N2O liquid [kmol]

    % Forward Difference Time Loop
    tf = 9.5;           % final time [s]
    tstep = 0.0005;     % time step [s]
    
    i_i = 0;
    i_f = tf/tstep;

    for i=i_i:i_f
        t = i*tstep;

        % Given functions of temperature:
        Vhat_l = Q2^(1+(1-To/Q3)^Q4)/Q1;
            %molar specific volume of liquid N2O [m^3/kmol]
        CVhat_He = C1 + C2*To + C3*To^2 + C4*To^3 + C5*To^4 - R;
            %specific heat of He at constant volume [J/(kmol*K)]
        CVhat_g = D1 + D2*((D3/To)/sinh(D3/To))^2 + D4*((D5/To)/cosh(D5/To))^2 - R;
            %specific heat of N2O gas at constant volume [J/(kmol*K)]
        CVhat_l = E1 + E2*To + E3*To^2 + E4*To^3 + E5*To^4;
            %specific heat of N2O liquid at constant volume, approx. same as at
            %constant pressure [J/(kmol*K)]
        Tr = To/Tc;  % reduced temperature
        delta_Hv = J1*(1 - Tr) ^ (J2 + J3*Tr + J4*Tr^2);    % heat of vapourization of N2O [J/kmol]
        P_sat = exp(G1 + G2/To + G3*log(To) + G4*To^G5);    % vapour pressure of N2O [Pa]
        dP_sat = (-G2/(To^2) + G3/To + G4*G5*To^(G5-1)) * exp(G1 + G2/To + G3*log(To) + G4*To^G5);
            %derivative of vapour pressure with respect to temperature
        Cp_T = (4.8 + 0.00322*To)*155.239;                  % specific heat of tank, Aluminum [J/(kg*K)]
        
        nitrous_density = 44.013/Vhat_l;
        
        % Calculations for Burnell Injector Model
        % Calculation of Burnell C-Coefficient
        C_Burnell = 0.315*exp(-0.00104*P_sat*0.000145038);
        % Calculation of choked flow nitrous mass flow rate
        mdot_Burnell = Cd*Ainj*sqrt(2*nitrous_density*C_Burnell*P_sat);
        molar_rate_Burnell = -mdot_Burnell/44.013;

        % Lots of equations from Margeret Mary's Masters Thesis
        % Used for calculation of pressure and temperature derivatives in
        % the nitrous tank
        % TODO: Go over this to make sure no other changes need to be made
        % to compensate for the burnell choked flow model
        % Simplified expression definitions for solution
        P = (n_He + n_go)*R*To/(V-n_lo*Vhat_l);
        a = m_T*Cp_T + n_He*CVhat_He + n_go*CVhat_g + n_lo*CVhat_l;
        b = P*Vhat_l;
        e = -delta_Hv + R*To;
        % f is the molar flow rate of the nitrous, if it is choked flow
        % (Chamber pressure is less than 0.8*nitrous vapor pressure), we
        % need to change f to be the choked flow value given by burnell
        f = -Cd*Ainj*sqrt(2/MW2)*sqrt((P-Pe)/Vhat_l);
        if(Pe < P*0.8)
            f = molar_rate_Burnell;
        end
        j = -Vhat_l*P_sat;
        k = (V - n_lo*Vhat_l)*dP_sat;
        m = R*To;
        q = R*n_go;

        Z = (-f*(-j*a+(q-k)*b))/(a*(m+j) + (q-k)*(e-b));
        W = (-Z*(m*a + (q-k)*e))/(-j*a + (q-k)*b);

        % Derivative Functions
        dT = (b*W+e*Z)/a;
        dn_g = Z;           %Rate (molar) of change of gaseous oxidizer [mol/s]
        if(Pe < P*0.8)
            dn2draindt = -molar_rate_Burnell;  %total molar drain rate of oxidizer [mol/s]
            dn_l = -dn_g - dn2draindt;           %Rate (molar) of change of liquid oxidizer [mol/s]
            doxdt = mdot_Burnell;    %total drain rate of oxidizer [kg/s]
        else
            dn_l = W;                   %Rate (molar) of change of liquid oxidizer [mol/s]
            dn2draindt = -dn_g - dn_l;  %total molar drain rate of oxidizer [mol/s]
            doxdt = dn2draindt * 44.013;    %total drain rate of oxidizer [kg/s]
        end
        
        %Oxidizer cannot return to the tank, should make this throw some
        %kind of error
        if(doxdt < 0) 
            doxdt = 0;
        end

        %regression rate of the fuel grain [m/s] rdot = aG^n, G =
        %doxdt/Aport
        drdt = a_coeff*realpow(doxdt/(pi*D*D/4),n_coeff);
        %Fuel mass flow rate [kg/s]
        dFdt = Fuel_Density * pi * D * Length * drdt;

        %O/F ratio
        of = doxdt/dFdt;

        %interpolated data
        C_Star = interp2(data_gridx,data_gridy,data_Cstar,of,Pe);
        Thrust_Coeff = interp2(data_gridx,data_gridy,data_thrust,of,Pe);
        Gamma = interp2(data_gridx,data_gridy,data_gamma,of,Pe);
        Temp = interp2(data_gridx,data_gridy,data_temp,of,Pe);
        %Current chamber volume (should calculate the constant to add to
        %the cylindrical section)
        Vc = 0.00011 + pi*D*D*Length/4;
        
        %Instantaneous change in chamber pressure
        dPe = (doxdt+dFdt) - Pe*A_Star*sqrt(Gamma/(R1*Temp))...
            *realpow(2/(Gamma+1),(Gamma+1)/(2*(Gamma - 1)));
        dPe = dPe * R1 * Temp / Vc;

        %Thrust [N]
        thrust = (doxdt + dFdt) * C_Star * Nozzle_Efficiency * Reaction_Efficiency * Thrust_Coeff;
        
        % Record variables for each time step in an array
        T(i+1,1) = t;
        T(i+1,2) = To;
        n_g(i+1,1) = t;
        n_g(i+1,2) = n_go;
        n_l(i+1,1) = t;
        n_l(i+1,2) = n_lo;
        Pres(i+1,1) = t;
        Pres(i+1,2) = P;
        PSat(i+1,1) = t;
        PSat(i+1,2) = P_sat;
        PE(i+1,1) = t;
        PE(i+1,2) = Pe;
        Dn_g(i+1,1) = t;
        Dn_g(i+1,2) = dn_g;
        Dn_l(i+1,1) = t;
        Dn_l(i+1,2) = dn_l;
        Dn_comb(i+1,1) = t;
        Dn_comb(i+1,2) = dn2draindt;
        Dn_ox(i+1,1) = t;
        Dn_ox(i+1,2) = doxdt;
        Dn_F(i+1,1) = t;
        Dn_F(i+1,2) = dFdt;
        Dn_m(i+1,1) = t;
        Dn_m(i+1,2) = doxdt + dFdt;
        Thrust(i+1,1) = t;
        Thrust(i+1,2) = thrust;
        OF(i+1,1) = t;
        OF(i+1,2) = of;
        Diam(i+1,1) = t;
        Diam(i+1,2) = D;
       
        %Forward Difference Method
        To = To + dT*tstep;
        n_go = n_go + dn_g*tstep;
        n_lo = n_lo + dn_l*tstep;
        Pe = Pe + dPe*tstep;
        D = D + 2*drdt*tstep;

        % Physical stops to kick out of loop
        if Pe>=P        %Chamber pressure cannot be greater than tank pressure (explosion)
            break
        end
        if n_lo<=0      % Stop when no liquid oxidizer remains
                        % We should also add a model for determining
                        % parameters when in gaseous flow
            break
        end
    end
    
    %Add in final values, send thrust and stuff to 0
    t = t + tstep;
    
    T(i+1,1) = t;
    T(i+1,2) = To;
    n_g(i+1,1) = t;
    n_g(i+1,2) = n_go;
    n_l(i+1,1) = t;
    n_l(i+1,2) = n_lo;
    Pres(i+1,1) = t;
    Pres(i+1,2) = P;
    PSat(i+1,1) = t;
    PSat(i+1,2) = P_sat;
    PE(i+1,1) = t;
    PE(i+1,2) = Pe;
    Dn_g(i+1,1) = t;
    Dn_g(i+1,2) = 0;
    Dn_l(i+1,1) = t;
    Dn_l(i+1,2) = 0;
    Dn_comb(i+1,1) = t;
    Dn_comb(i+1,2) = 0;
    Dn_ox(i+1,1) = t;
    Dn_ox(i+1,2) = 0;
    Dn_F(i+1,1) = t;
    Dn_F(i+1,2) = 0;
    Dn_m(i+1,1) = t;
    Dn_m(i+1,2) = 0;
    Thrust(i+1,1) = t;
    Thrust(i+1,2) = 0;
    OF(i+1,1) = t;
    OF(i+1,2) = 0;
    Diam(i+1,1) = t;
    Diam(i+1,2) = D;
    
    % Plot thrust curve (should change Newtons to lbf)
    figure(11), plot(Thrust(:,1),Thrust(:,2), 'r','LineWidth',1),grid, ...
            title('Thrust vs. Time'), ...
            xlabel('Time [s]'), ...
            ylabel('Thrust [N]'),drawnow;

    % Calculate total impulse using trapezoidal approximation of integral
    imp = trapz(Thrust(:,1),Thrust(:,2));
    
    % Output various other parameters (Should change to imperial units)
    global OUTPUT_RUN_GRAPHS;
    
    if OUTPUT_RUN_GRAPHS == 1
        figure(6), plot(T(:,1),T(:,2), 'r','LineWidth',2),grid, ...
            title('Temperature vs. Time'), ...
            xlabel('Time [s]'), ...
            ylabel('Temperature [K]'),drawnow;
        figure(7), plot(n_g(:,1),n_g(:,2)*44.013*2.20462,'b',n_l(:,1),n_l(:,2)*44.013*2.20462,'g','LineWidth',2),grid, ...
            title('N2O vs. Time'), ...
            xlabel('Time [s]'), ...
            ylabel('lb of N2O [lb]'), ...
            legend('lb of N2O gas','lb of N2O liquid'),drawnow;
        figure(8), plot(Pres(:,1),Pres(:,2)*0.000145038,'m',...
            PE(:,1),PE(:,2)*0.000145038,'c',PSat(:,1),PSat(:,2)*0.000145038,'b','LineWidth',2),grid, ...
            title('Pressure vs. Time'), ...
            xlabel('Time [s]'), ...
            ylabel('Pressure [Psi]'), ...
            legend('tank pressure','chamber pressure'),drawnow; 
        figure(9), plot(Dn_g(:,1),Dn_g(:,2),'m',Dn_l(:,1),Dn_l(:,2),'c',Dn_comb(:,1),Dn_comb(:,2), 'r','LineWidth',2),grid, ...
            title('Molar Flow Rate vs. Time'), ...
            xlabel('Time [s]'), ...
            ylabel('Rate [kmol/s]'), ...
            legend('Gas','Liquid','Combined'),drawnow;
        figure(10), plot(Dn_ox(:,1),Dn_ox(:,2),'m',Dn_F(:,1),Dn_F(:,2),'c',Dn_m(:,1),Dn_m(:,2), 'r','LineWidth',2),grid, ...
            title('Mass Flow Rate vs. Time'), ...
            xlabel('Time [s]'), ...
            ylabel('Rate [kg/s]'), ...
            legend('Oxidizer','Fuel','Combined'),drawnow;
        figure(11), plot(Thrust(:,1),Thrust(:,2)*0.224809, 'r','LineWidth',1),grid, ...
            title('Thrust vs. Time'), ...
            xlabel('Time [s]'), ...
            ylabel('Thrust [lbf]'),drawnow;
        figure(12), plot(OF(:,1),OF(:,2), 'r','LineWidth',1),grid, ...
            title('OF vs. Time'), ...
            xlabel('Time [s]'), ...
            ylabel('O/F'),drawnow;
        figure(13), plot(Diam(:,1),Diam(:,2), 'r','LineWidth',1),grid, ...
            title('Diameter vs. Time'), ...
            xlabel('Time [s]'), ...
            ylabel('Diameter [m]'),drawnow;
    end
end

%make sure that we don't put too much N2O into the tank (unused atm, should fix and use)
function good = checkICs(m_loaded,Ti,V)
    good = true;
    % Perry's Chemical Engineers' Handbook Property Equations
    G1 = 96.512;        % vapour pressure of N2O [Pa] coefficients
    G2 = -4045;         % valid for Temp range [182.3 K - 309.57 K]
    G3 = -12.277;      
    G4 = 2.886e-5;
    G5 = 2;

    Q1 = 2.781;         %molar specific volume of liquid N2O [m^3/kmol] coefficients
    Q2 = 0.27244;
    Q3 = 309.57;
    Q4 = 0.2882;

    % Initial conditions
    n_to = m_loaded/MW2;                                        % initial total N2O in tank [kmol]
    Vhat_li = Q2^(1+(1-Ti/Q3)^Q4)/Q1;                           % molar volume of liquid N2O [m^3/kmol]    
    To = Ti;                                                    % initial temperature [K]   
    P_sato = exp(G1 + G2/To + G3*log(To) + G4*To^G5);           % initial vapour pressure of N2O [Pa]
    n_go = P_sato*(V - Vhat_li*n_to)/(-P_sato*Vhat_li + R*To);  % initial N2O gas [kmol]
    n_lo = (n_to*R*To - P_sato*V)/(-P_sato*Vhat_li + R*To);     % initial N2O liquid [kmol]
    
    if n_go < (n_to/52.0)
        good = false;
    end
end