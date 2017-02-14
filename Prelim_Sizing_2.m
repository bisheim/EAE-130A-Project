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
TSFC_initial = (1+0.35*Mach)*sqrt(theta) % Low Power Turbofan Settings: Mattingly
V = Mach * a_std * sqrt(theta); % Velocity in ft/s
range = 4000 * 6080; % Required flight range in ft
Loiter = 2; % hours

W_crew = 200;
Num_crew = 2; % Fixed
W_crew = W_crew*Num_crew;

W_1pass = 230; % Weight of one pass
Num_pass = 19; % VARIABLE
W_payload = (W_1pass)*Num_pass; 

W_tfo = 0; % Assume zero for now
W_fuelmax = 80000;
W_empty = 50000; % Aerion

%% Sadreay
W_fuel = 60000; % INitial Estimate
W_to = W_fuel + W_crew + W_empty + W_payload %Initial Estimate


    W_to_1 = (W_payload+W_crew)/(1-(W_fuel)-(W_empty)); % Fraction

        RF = L_overD*V*3600/TSFC_initial; % Range Factor
        c_ratio = exp(-range/RF); % Cruise Range
        l_ratio = exp(-TSFC_initial*Loiter/L_overD); % Loiter ratio 
    W6_1 = 0.98*0.97*0.99*0.997*c_ratio * l_ratio;

    W_fuel = (1-(W6_1))*1.05 % Fuel fraction Wof fuel div by wto 
    W_empty = (1.13*10^(-6)*W_to + 0.48);
    
W_to = W_fuel + W_crew + W_empty + W_payload 








