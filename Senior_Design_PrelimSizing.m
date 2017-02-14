%% Senior Design Preliminary Sizing
% EAE 130A 
% 2.8.16
% 
clear all
close all
clc
% Assign Values
%% Weight of people
W_1pass = 180; % Weight of one pass
W_pluggage = 50;
W_crew = 200;
Num_pass = 19; % VARIABLE
Num_crew = 2; % Fixed

W_payload = Num_pass*(W_1pass+W_pluggage)+W_crew*Num_crew;

%%
L_overD = 8; % FROM SADRAEY

Altitude = 50000; % Below top of topopause
Mach = 1.6; % Minimum Cruise Requirement
theta = 0.7519; % DETERMINED FROM ALTITUDE Tables
a_std = 1116; % Speed of sound 
TSFC_initial = (1+0.35*Mach)*sqrt(theta); % Low Power Turbofan Settings: Mattingly
V = Mach * a_std * sqrt(theta); % Velocity in ft/s
range = 4000 * 6080; % Required flight range in ft
Loiter = 2; % hours


%% Nikos Data Get initial value for take off weight W_to
W_crew = W_crew*Num_crew;
W_empty = 50000; % Aerion
W_1pass = 180; % Weight of one pass
W_pluggage = 50;
Num_pass = 19;
W_payload = (W_1pass + W_pluggage)*Num_pass; 
W_tfo = 0; % Assume zero for now
W_fuelmax = 80000;
%W_reserve = 1,000; % TOTALLY JUST PULLED THIS OUT OF THE AIR! JUSTIFY

W_oe = W_empty + W_crew+W_tfo;

W_to = W_oe + W_fuelmax + W_payload
%W_to

    %% Sadreay
W_fuel = 60000;
W_to = (W_payload+W_crew)/(1-(W_fuel/W_to)-(W_empty/W_to))
 
RF = L_overD*V*3600/TSFC_initial; % Range Factor
c_ratio = exp(-range/RF); % Cruise Range
l_ratio = exp(-TSFC_initial*Loiter/L_overD);
W6 = W_to*0.98*0.97*0.99*0.997*c_ratio * l_ratio

W_fuel = W_to * (1-(W6/W_to))*1.05; % Including extra fuel
W_empty = W_to*(1.13*10^(-6)*W_to + 0.48)

    
%%















%% Use initial value to recalculate the taekoff weight
W_fuel = 60000;
itermax = 100;
iter = 1;
resmin = 1;
res = 100;

for i=1:100
            %W_ref = W_fuel;
            %W_oe = W_empty + W_crew+W_tfo;
            %W_to = W_oe + W_fuel + W_payload;
            W1 = W_to;
            W2 = 0.97*W1;% Taxi to Climbout- Empirical value from class 
            W3 = 0.97*W2; % Accelerate and climb to cruise
            RF = L_overD*V*3600/TSFC_initial; % Range Factor
            W4 = W3 * exp(-range/RF); % Cruise Range
            W5 = W4 /(exp(TSFC_initial*Loiter/L_overD)); % Loiter Equation for loiter at airport 
            W6 = 0.99*W5; % Descent
            W7 = 0.997*W6; % Approach and Landing
            W_fuel = (1-(W7/W1))*W_to + W_reserve;
            
            W_oe = W_empty + W_crew+W_tfo;
            W_to = W_oe + W_fuel + W_payload;
            W_empty = W_to*(1.13*10^(-6)*W_to + 0.48);

            plot(i,W_fuel,'o')
            %plot(i,W4,'r+')
            title('Weight Convergence')
            xlabel('Iteration'); ylabel('Weight of Fuel (lbs)');
            hold on
            %res = W_ref - W_fuel;
            iter = iter+1;
end    
    
    %% Sadreay

W_to = (W_payload+W_crew)/(1-(W_fuel/W_to)-(W_empty/W_to))
 
RF = L_overD*V*3600/TSFC_initial; % Range Factor
W_cratio = exp(-range/RF); % Cruise Range
W6 = W1*0.98*0.97*0.99*0.997*W_cratio * exp(-TSFC_initial*Loiter/L_overD);

W_fuel = W_to * (1-(W6/W_to))*1.05; % Including extra fuel
W_empty = W_to*(1.13*10^(-6)*W_to + 0.48)

    
    