%% Preliminary Sizing Part 
% 2/13/17
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
Loiter = 0.5; % hours

W_crew = 200;
Num_crew = 2; % Fixed
W_crew = W_crew*Num_crew;

W_1pass = 230; % Weight of one pass
Num_pass = 19; % VARIABLE
W_payload = (W_1pass)*Num_pass; 

W_tfo = 0; % Assume zero for now
W_fuelmax = 80000;
W_empty = 50000; % Aerion
W_to = 100000; % Aerion
RF = L_overD*V*3600/TSFC_initial; % Range Factor


%% Running it for Different Numbers of Passengers
for j = 6:19
W_payload = (W_1pass)* j; 

        for i=1:50
        W1 = W_to;
        W2 = 0.97*W1;% Taxi to Climbout- Empirical value from class 
        W3 = 0.97*W2; % Accelerate and climb to cruise
                        RF = L_overD*V*3600/TSFC_initial; % Range Factor
        W4 = W3 * exp(-range/RF); % Cruise Range
        W5 = W4 * (exp(-TSFC_initial*Loiter/L_overD)); % Loiter Equation for loiter at airport 
        W6 = 0.99*W5; % Descent
        W7 = 0.997*W6; % Approach and Landing
        
        %W_empty = (1.13 * 10^(-6) * W_to^2) + (0.48*W_to)
        W_fuel = (1-(W7/W1))*W_to*1.05;
        W_to = W_fuel + W_payload + W_empty+ W_crew +W_tfo;
        %figure(j-5)
        %plot(i,W_to, 'ro'); hold on; title('Weight Sizing'); xlabel('iteration');ylabel('Take Off weight in lbs')

        flight_prof=[W1,W2,W3,W4,W5,W6,W7];
        end

    weight_takeoff(j-5) = W_to;
    weight_fuel(j-5) = W_fuel;

end
x=[6:19];
figure(20); hold on
plot(x,weight_takeoff,'ro')
plot(x,weight_fuel,'b+')
legend('Total Takeoff Weight', 'Fuel Weight','Location','Best');
xlabel('Number of Passengers'); ylabel('Weight in lbs')



