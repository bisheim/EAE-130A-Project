%% Senior Design Preliminary Sizing
% EAE 130A 
% 2.8.16
% 
clear all
close all
clc

%%
L_overD = 8; % FROM SADRAEY
Altitude = 50000; % Below top of topopause
Mach = 1.6; % Minimum Cruise Requirement
theta = 0.7519; % DETERMINED FROM ALTITUDE Tables
a_std = 1116; % Speed of sound 
TSFC_initial = (1+0.35*Mach)*sqrt(theta); % Low Power Turbofan Settings: Mattingly
V = Mach * a_std * sqrt(theta); % Velocity in ft/s
range = 4000 * 6080; % Required flight range in ft
Loiter = .75; % hours

W_crew = 200;
Num_crew = 2; % Fixed
W_crew = W_crew*Num_crew;

W_1pass = 230; % Weight of one pass
Num_pass = 19; % VARIABLE
W_payload = (W_1pass)*Num_pass; 

W_tfo = 0; % Assume zero for now
W_fuelmax = 80000;
% W_empty = 50000; % Aerion

%% Sadreay
% W_fuel = 60000; % INitial Estimate
W_to = 150000; %guess Initial Estimate

RF = L_overD*V*3600*TSFC_initial; % Range Factor
c_ratio = exp(-range*.8/3600/V/L_overD); % Cruise Range Use .8 for TSFC according to Sandray
l_ratio = exp(-.7*Loiter/L_overD); % Loiter ratio  Use .7 for TSFC according to Sandray

%empirical from Sandray
W_takeoff = .98; 
W_climb = 0.97;
W_descent = 0.99; 
W_approachland = 0.997;

W6_1 = W_takeoff*W_climb*W_descent*W_approachland * c_ratio * l_ratio;
W_fuel = (1-(W6_1))*1.06; % Fuel fraction Wof fuel div by wto 

error = 15; % percent error initlal sizing
%%
while abs(error) > .0001
    W_to_temp = W_to;
    W_empty = 1.02*(W_to_temp)^(-.06); %estimate from Aircraft design a Conceptual Approach (1992) by DANIEL P RAYMER
    W_to = (W_payload+W_crew)/(1-(W_fuel)-(W_empty)); % Fraction
    
    error = (W_to-W_to_temp)/W_to_temp*100;
end
%disp(W_to)
%disp(error)
%% Find Wing Area and Thrust
%Stall Speed
rho_o = 0.002378; %slug/ft^3 @sea level
V_landing = 160*1.68781; %ft/s, Concorde Landing Speed
V_stall = V_landing/1.3;
C_Lmax = 2.0; %Design Parameter (Research this Weekend)
W_Stall = 0.5*rho_o*V_stall^2*C_Lmax; %lbf/ft^2

%Max Speed
W_S_var = 1:200;
C_do = .04; %Design Parameter (Research this Weekend)
a = rho_o*C_do/2;
sigma_cruise = 0.1151/0.7519; %Mattingly Tables
rho_cruise = sigma_cruise*0.0767/32.174; %slug/ft^3
e = 0.7; %Design Parameter
AR = 2.7563; %Design Parameter (Research)
K = 1/(pi*e*AR);
b = 2*K/(rho_cruise*sigma_cruise);
V_cruise = V;
V_max = 1.2*V_cruise;
T_W_maxspeed = a*V_max^2./W_S_var + b/V_max^2*(W_S_var);

%Takeoff Run Requirements : Saadraey Page 
mu = 0.05; % Coeff Friction runway
g = 32.174; % Gravity
S_to = 7000; %Runway length from RFP
C_Dlg_o = 0.012; % Landing Gear dra coeff S 4.69b
C_Dhild_to = 0.003; % Hi lift devices on TO withou Kreuger flaps S. 4.96b
C_Lc =0.05; % Aircraft Cruise Lift Coeff
C_Lflap_to = 0.5; % Lift coeff from flaps at takeoff

V_r = 1.1*V_stall; % Rotational speed at wheels up
%C_LR = ((2*g)/(rho_o*V_r^2)).*(W_S_var); %Rotational Lift Coefficient
C_LR = C_Lmax/(1.15^2); % Rotational Lift Coefficient
C_Lto = C_Lc + C_Lflap_to;
C_Dto_o = C_do + C_Dlg_o+C_Dhild_to;
C_Dto = C_Dto_o+K*C_Lto^2; 
C_DG = C_Dto-mu*C_Lto;

T_W_takeoffrun = (mu-(mu+(C_DG./C_LR)).*(exp(0.6*rho_o*g*C_DG*S_to.*(W_S_var).^(-1))))./(1-exp(0.6*rho_o*g*C_DG*S_to.*(W_S_var).^(-1)));

% Rate of Climb Requirements
ROC = 6000; %  Design Parameter (Research )fpm Concorde was 5000 
T_W_ROC = (ROC/sqrt((2/(rho_o*sqrt(C_do/K)))))+1/(L_overD);

%Ceiling Requirements
ROC_cruise = 300; %fpm cruise ceiling 
T_W_cruise_ceiling = ROC_cruise./(sigma_cruise*sqrt(2./(rho_cruise*sqrt(C_do/K)).*W_S_var))+(1/(sigma_cruise*L_overD));


plot(W_S_var,T_W_maxspeed,'r'); hold on
plot(W_S_var,T_W_takeoffrun,'g');
plot(W_S_var,T_W_cruise_ceiling,'k');
plot(W_Stall,1:160,'b');
plot(W_S_var,T_W_ROC,'c');
legend('Max Speed','Takeoff', 'Ceiling', 'Stall', 'ROC', 'Location', 'best')
