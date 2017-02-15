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
W_to
error