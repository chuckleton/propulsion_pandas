%Henrique Bittencourt Netto Monteiro - hbm2122
%Columbia Space Initiative, New York, NY
%Fall 2018

% CODE TO CALCULATE:

   % THERMOFLUID PARAMETERS
      % Ratio of throat/exit velocity
      % Ratio of exit/throat area
      % Throat/Chamber pressure ratio
      % Critical pressure at throat
      % Sonic velocity at throat
      % Mass flow at nozzle
   
   % MODELLING PARAMETERS  
       % Nozzle equation


%calculating the ratio between the throttle area and the exit area
disp(" ")
disp("k = 1.1624     P_exit = 1atm / 14.696psi     P_chamber = 325psi")
disp("R = 309.7 J/kg*K     T_inlet = 3211K");
disp(" ");
answer = input('Do you want to use the values above? Answer "y" or "n" ', 's');

if answer == 'y'
    k = 1.1624;
    pe = 14.696;
    pi = 325;
    R = 309.7;
    Tinl = 3211.2737;
else
    k = input('What is the specific heat ratio k? ');
    pe = input('What is the exit pressure? ');
    pi = input('What is the chamber pressure? ');
    Tinl = input('What is the temperature at the nozzle inlet? ');
    R = input('What is the gas constant? ');
end

%calculating terms for area ratio
term1 = ((k+1)/2)^(1/(k-1));
term2 = (pe/pi)^(1/k);
term3 = (k+1)/(k-1);
term4 = 1 - (pe/pi)^((k-1)/k);

%calculating throat/exit velocity ratio
vel_ratio = sqrt(term3*term4);
%calcuting exit/throat area ratio
area_ratio = 1/(term1*term2*vel_ratio);
disp(" ")
disp("Ratio of throat/exit velocity: " + num2str(vel_ratio));
disp("Ratio of exit/throat area: " + num2str(area_ratio));

%calculating pressure at throat for the flow to be Mach 1
p_ratio = [2/(k+1)]^(k/(k-1));
%calcuting velocity at throat for critical pressure
thr_vel = sqrt((2*k/(k+1))*R*Tinl);
%calculating the mass flow
term5 = (2/(k+1))^((k+1)/(k-1));
term6 = sqrt(k*R*Tinl);
pi_SI = pi*6894.757;
m_flowXX = pi_SI*k*(sqrt(term5)/term6);
disp("Throat/Chamber pressure ratio: " + num2str(p_ratio));
disp(" ")
disp("Critical pressure at throat: " + num2str(pi*p_ratio) + " psi");
disp("Sonic velocity at throat: " + num2str(thr_vel) + " m/s");
disp("Mass flow at nozzle: " + num2str(m_flowXX) + " * Area_throat kg/m2*s");
disp(" ")

length = input('Length of nozzle: ');
dim2 = input('Diameter of exit: ');

dim1 = dim2/sqrt(area_ratio);
A = (dim2/2 - dim1/2)/(length^(1/2));
disp("Nozzle equation: " + num2str(dim1/2) + " + " + num2str(A) + "*x^(1/2)");
disp(" ")


