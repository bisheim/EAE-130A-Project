%% Senior Design Preliminary Sizing
% EAE 130A 
% 2.8.16
% 
clear all
close all
clc

%%
L_overD = 8.4; % FROM SADRAEY (Maureen)
Altitude = 50000; % Below top of topopause
Mach = 1.6; % Minimum Cruise Requirement
theta = 0.7519; % DETERMINED FROM ALTITUDE Tables
a_std = 1116; % Speed of sound 
TSFC_initial = (1+0.35*Mach)*sqrt(theta); % Low Power Turbofan Settings: Mattingly
V = Mach * a_std * sqrt(theta); % Velocity in ft/s
range = 4000 * 6080; % Required flight range in ft
Loiter = 1; % hours

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
W_S_var = 1:100;
e = 0.756; %Design Parameter (Nikos)
AR = 2.5; %Design Parameter (Nikos)
C_do = .04; %Design Parameter (Bond)
K = 1/(pi*e*AR);
C_do = .04; %Design Parameter (Bond)
rho_o = 0.002378; %slug/ft^3 @sea level
L_overD_max = 9.699;
a = rho_o*C_do/2;
sigma_cruise = 0.1151/0.7519; %Mattingly Tables
rho_cruise = sigma_cruise*0.076647/32.174; %slug/ft^3
V_cruise = V;


%Stall Speed at landing
rho_o = 0.002378; %slug/ft^3 @sea level
sigma_o = 1; % at 3000 ft altitude of the airport
V_landing = 120*1.68781; %ft/s, Spike S-512 Landing Speed
V_stall = 120*1.68781; %ft/s based on stall speed of Spike S-152
C_Lmax = 1.6;
W_Stall = 0.5*rho_o*sigma_o*V_stall^2*C_Lmax; %lbf/ft^2
W_Stall = [W_Stall W_Stall];
T_W_Stall = [0 100];
%flaps up
C_Lland = 1.8;
W_Stall_land = 0.5*rho_o*sigma_o*V_stall^2*C_Lland; %lbf/ft^2
W_Stall_land = [W_Stall_land W_Stall_land];
T_W_Stall_land = [0 100];

%Takeoff Run Requirements : Saadraey Page 131
s_tofl=7000;
T_W_takeoffrun = 37.5*W_S_var/s_tofl/sigma_o/(C_Lmax+0.2);

% Drag polars
a = -2.5229;% based on business jet correlation in Roskam
b = 1;% based on business jet correlation in Roskam
c = -1.1868; % based on supersonic aircraft Roskam
d = 0.9609;% based on supersonic aircraft Roskam
f = 10^(a+b*(c+d*log10(W_to)));

C_L_var = C_Lmax;
CD_clean = C_do+K*C_L_var.^2;% clean
CD_to_gu = (C_do+.01)+1/(pi*AR*e)*C_L_var.^2;% takeoff landing gear up
CD_to_gd = (C_do+.01+.025)+1/(pi*AR*e)*C_L_var.^2;% takeoff landing gear down
CD_l_gu = (C_do+.055)+1/(pi*AR*e)*C_L_var.^2;% landing landing gear up
CD_l_gd = (C_do+.055+.025)+1/(pi*AR*e)*C_L_var.^2;% landing landing gear down

CD_to_gu = (C_do)+1/(pi*AR*e)*C_L_var.^2;% takeoff landing gear up
CD_to_gd = (C_do)+1/(pi*AR*e)*C_L_var.^2;% takeoff landing gear down
CD_l_gu = (C_do)+1/(pi*AR*e)*C_L_var.^2;% landing landing gear up
CD_l_gd = (C_do)+1/(pi*AR*e)*C_L_var.^2;% landing landing gear down

%Far part 25 requirements, must also do 50F higher than standard day
CGR_OEI_to_gu = 0.012; %25.111
CGR_OEI_to_gd = 0;%25.121
CGR_OEI_to_gu_noge = 0.024; %no ground effect %25.121
CGR_OEI_climb_gu = 0.012;%25.121

CGR_AEI_l_gd = 0.032;%25.119
CGR_AEI_l_gd_approach = 0.021;%25.121

T_W_OEI_to_gu = 2*((L_overD)^-1+CGR_OEI_to_gu); % accounts for 50F with the /.8 parameter
T_W_OEI_to_gd = 2*((L_overD)^-1+CGR_OEI_to_gd); % accounts for 50F with the /.8 parameter
T_W_OEI_to_gu_noge = 2*((L_overD)^-1+CGR_OEI_to_gu_noge);
T_W_AEI_l_gd = 2*((L_overD_max)^-1+CGR_AEI_l_gd);
T_W_AEI_l_gd_approach = 2*((L_overD_max)^-1+CGR_AEI_l_gd_approach);

% Cruise
q = .5*rho_cruise*V_cruise^2;
C_D0_withcomp = C_do;%compressibility effects
T_W_cruise = (C_D0_withcomp*q./W_S_var+W_S_var/(q*e*AR));

% manuevering
q = .5*rho_cruise*V_cruise^2;
C_D0_withcomp = C_do;%compressibility effects
psi = 1.5; %standard turn rate of 1.5°/s CITE: https://www.paramountbusinessjets.com/aviation-terminology/standard-rate-turn.html
V_psi = 250*1.68781; %minimum velocity for turning at 1.5°/s
n_max = (V_psi*psi/32.174+1)^.5;
T_W_turn = C_D0_withcomp*q./W_S_var + W_S_var*n_max^2/(q*e*AR);

% Find the design point
index_dp  = sum(W_Stall(1)>W_S_var);
dp = (W_Stall(1)-W_S_var(index_dp))*(T_W_cruise(index_dp+1)-T_W_cruise(index_dp))/(W_S_var(index_dp+1)-W_S_var(index_dp))+T_W_cruise(index_dp);



%plotting
plot(W_S_var,T_W_cruise,'r'); hold on
plot(W_Stall,T_W_Stall, 'b');
plot(W_S_var,T_W_takeoffrun,'g');
plot([0 100],[T_W_OEI_to_gu_noge T_W_OEI_to_gu_noge],'k-');
plot(W_Stall_land,T_W_Stall_land, 'c');
plot([0 100],[T_W_AEI_l_gd T_W_AEI_l_gd],'m-');
plot(W_Stall(1), dp, 'ro', 'MarkerFaceColor', 'b','MarkerSize', 8)


ylim([0,1])
xlabel('W / S_{TO} (lbf/ft^2)')
ylabel('T / W_{TO}')
legend('Cruise','Stall', 'Takeoff Run', 'ROC OEI','Landing Stall','ROC AEI','Design Point', 'Location', 'northeast')

[c,h] = contourf(xtest,ytest,[0 1],[2 2]); 
hp = findobj(h,'type','patch'); 
hatchfill(hp);

%% Drag Polar

%CD = CD0 + CL^2/(pi AR e)

CL_dragpolar = 0:.01:0.5;
CD_dragpolar_AR25 = C_do + CL_dragpolar.^2/(pi*e*2.5);
CD_dragpolar_AR20 = C_do + CL_dragpolar.^2/(pi*e*2.0);
CD_dragpolar_AR30 = C_do + CL_dragpolar.^2/(pi*e*3.0);
CD_dragpolar_AR40 = C_do + CL_dragpolar.^2/(pi*e*4.0);
CD_dragpolar_AR50 = C_do + CL_dragpolar.^2/(pi*e*5.0);
CD_dragpolar_AR60 = C_do + CL_dragpolar.^2/(pi*e*6.0);
CD_dragpolar_AR80 = C_do + CL_dragpolar.^2/(pi*e*8.0);

CL_atcruise = W_Stall(1)*W_climb* c_ratio/(.5*rho_cruise*V_cruise^2);

%plotting
figure()
plot(CD_dragpolar_AR20,CL_dragpolar, 'm'); hold on
plot(CD_dragpolar_AR25,CL_dragpolar, 'c');
plot(CD_dragpolar_AR30,CL_dragpolar, 'r');
plot(CD_dragpolar_AR40,CL_dragpolar, 'g');
plot(CD_dragpolar_AR50,CL_dragpolar, 'b');
plot(CD_dragpolar_AR60,CL_dragpolar, 'k');
plot(CD_dragpolar_AR80,CL_dragpolar, '--k');


xlabel('C_D')
xlim([0.036,0.05])
ylabel('C_L', 'rotation',90)
ylim([0,0.3])
set(get(gca,'YLabel'),'Rotation',0)
legend('AR = 2.0','AR = 2.5', 'AR = 3.0', 'AR = 4.0',...
    'AR = 5.0','AR = 6.0','AR = 8.0', 'Location', 'best')

%%

% %Max Speed
% C_do = .04; %Design Parameter (Bond)
% a = rho_o*C_do/2;
% sigma_cruise = 0.1151/0.7519; %Mattingly Tables
% rho_cruise = sigma_cruise*0.076647/32.174; %slug/ft^3
% e = 0.756; %Design Parameter (Nikos)
% AR = 2.5; %Design Parameter (Nikos)
% K = 1/(pi*e*AR);
% b = 2*K/(rho_cruise*sigma_cruise);
% V_cruise = V;
% V_max = 1.2*V_cruise;
% T_W_maxspeed = a*V_max^2./W_S_var + b/V_max^2*(W_S_var);
% 
% 
% 
% %climb
% ne=2;
% AEO=ne/(ne-1);
% L_overD_max = 9.699;
% T_W_climb = (AEO)*(1/L_overD_max+0.03);
% W_climb = [0 100];
% T_W_climb = [T_W_climb T_W_climb];
% 
% % Rate of Climb Requirements
% ROC = 5000; %  Design Parameter (Research )fpm Concorde was 5000 
% L_overD_max = 9.699; % Design parameter (Maureen)
% T_W_ROC = ((ROC./sqrt((2/(rho_cruise*sqrt(C_do/K))).*W_S_var))+1/(L_overD_max))*.1151/.7519;
% 
% %Ceiling Requirements
% ROC_service = 100; %fpm cruise ceiling 
% sigma_service = 0.2967*exp(1.7355-4.8075*10^-5*55000);
% T_W_service_ceiling = ROC_service./(sigma_service*sqrt(2./(rho_o*sigma_service*sqrt(C_do/K)).*W_S_var))+(1/(sigma_service*L_overD_max));
% 
% %Wing loading for cruise from Raymer pg 93
% 
% 
% 
% 
% % plot
% plot(W_S_var,T_W_maxspeed,'r'); hold on
% plot(W_S_var,T_W_takeoffrun,'g');
% plot(W_S_var,T_W_service_ceiling,'k');
% plot(W_S_var,T_W_ROC,'m');
% plot(W_Stall,T_W_Stall, 'b');
% plot(W_climb,T_W_climb);
% ylim([0,8])
% legend('Max Speed','Takeoff', 'Ceiling', 'ROC', 'Stall', 'Location', 'best')
