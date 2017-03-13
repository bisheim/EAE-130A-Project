%% Preliminary weight and Balance
clear all; close all; clc;




% Find average ratios of component weight fractions for a GW. 
GW = 122580; % lb Take off weight
structure = [0.242 0.310 .310 .284];
%structure_ave = mean(structure);
powerplant = [0.074 0.129 0.085];
%powerplant_ave = mean(powerplant);
fixed_equip = [0.090 0.119 0.164];
empty_weight = [0.406 0.562 0.495];
wing = [0.095 0.128 0.103 0.110];
empennage = [0.020 0.023 0.029 0.028];
h_tail = [0.015189];
v_tail = [0.007468];
fuselage = [0.073 0.093 0.122 0.108];
nacelle = [0.017 0.016 0.015];
landing_gear = [0.037 0.051 0.040 0.038];

ave_sizing = [mean(structure) mean(powerplant) mean(fixed_equip) mean(empty_weight) mean(wing) mean(empennage) mean(fuselage) mean(nacelle) mean(landing_gear)];
weight_val = ave_sizing .* GW;

%Weight_frac = [mean(fuselage) mean(wing) mean(h_tail) mean(v_tail) mean(empennage) mean(powerplant) mean(landing_gear) mean(fixed_equip) 0.0028 0.0033 0.456 0.036];
Weight_frac = [mean(fuselage) mean(wing) mean(h_tail) mean(v_tail) mean(powerplant) mean(landing_gear) mean(fixed_equip) 0.0028 0.0033 0.456 0.036];

Comp_weight = Weight_frac .* GW;

%% Find the overall locations of each component
    % Start with reference locations
A_x = 0 +20; % Nose
A_z = 0 +20; % Z nose
B_x = 70 +20; % front tip of wing
B_z = -0.80 +20; % z wing
C_x = 125 +20; % front of vert stab
C_z = 2.75 +20;% z vertstab
D_x = 138.0 +20; % engine inlet
D_z = 4.0 +20; % z engine
E_x = 29.643 +20; % front of cockpit
E_z = -2.8 +20; % Z of cockpit
F_x = 47.5 +20; % front of cabin
F_z =  -2.00 +20; %z cabin
G_x = 25.504 +20; % front tip of canard
G_z = 3 +20; % z canard

%% Componenet cg locations
    % all starred quantities must be fillled in from geometry
% Geometry that needs to be defined
length_fuselage = 150;%**** in feet for twin engine aft mounted fuselage

c_root = 52.4720;%*** Root chord wing
le_sweep = 60; % *** Leading edge sweep angle wing
te_sweep = 19.10284; % *** Trailing edge sweep wing
b = 75.73688;% **** Wing span b    


c_root_can = 25.95758;%*** Root chord wing cannard
le_sweep_can = 60; % *** Leading edge sweep angle wing cannard
te_sweep_can = 19.10284; % *** Trailing edge sweep wing cannard
b_can = 4.00;% **** cannard span b


c_root_tail = 15;%*** Root chord of vertical stabilizer
vle_sweep = 40; % *** Leading edge sweep angle of vertical stablizer
vte_sweep = 22.90080; % *** Trailing edge sweep of vertical stabilizer
v_b = 12;% **** vertical tail length

length_engine = 12.8300; % Length of engine from P&W specs
diam_engine = 4.1; %Engine diameter in front
length_cockpit = 10; %****
diam_front = 7.2; % **** diameter of fuselage at front 
diam_rear = 9; % Diameter of fuselage at rear
length_cabin = 54;% *****

%% Calculate x cg locations
% Fueslage: reference point, nose of aircraft
x_fuselage = 0.485 * length_fuselage;
x_cg_fuselage = x_fuselage + A_x;

z_cg_fuselage = A_z;

% wing: reference point is the leading tip of the swept wing
x_wing = 0.70*(1-0.35)*(b/2)*(tand(le_sweep)-tand(te_sweep)) + (tand(le_sweep)*(0.35*(b/2)));
x_cg_wing = x_wing + B_x;

z_cg_wing = B_z;

% Empennage: reference point: leading edge of tail root 
x_vtail = 0.42*(1-0.38)*v_b*(tand(vle_sweep)-tand(vte_sweep)) + (tand(vle_sweep)*(0.38*(v_b))); % for no horizontal sabilizer
%x_cg_vert_stab = 0.42*(1-0.55)*v_b*(tand(vle_sweep)-tand(vte_sweep)) +
    %(tand(vle_sweep)*(0.55*(v_b))) % if T tail configuration
x_cg_vtail = x_vtail + C_x;

z_cg_vtail = x_vtail + C_z;

% Canard: reference point is the leading tip of the swept wing
x_canard = 0.70*(1-0.35)*(b_can/2)*(tand(le_sweep_can)-tand(te_sweep_can)) + (tand(le_sweep_can)*(0.35*(b_can/2)));
x_cg_canard = x_canard + G_x;

z_cg_canard = G_z;

%Engine group : reference point inlet of engine **** needs revision
x_engine = 0.5*length_engine;
x_cg_engine = x_engine + D_x;

z_cg_engine = D_z;

%Landing gear: reference point, nose of the aircraft
    %assume mostly vertical struts so cg is where the landing gear are
    %placed
x_landgear = 76; %****
x_cg_landgear = x_landgear + A_x;

z_cg_landgear = -5 +20;

% Fixed Equipment: ref point nose of aircraft
x_fixed_equp = 0.3 * length_fuselage; % *** Needs Far Mosre Detail
x_cg_fixed_equp = x_fixed_equp + A_x;

z_cg_fixed_equp = A_z;

% Trapped fuel and oil: ref point, nose of aircraft
    % near engine case
x_tfo = 0; %**** assuming its near engine inlet
x_cg_tfo = x_tfo + D_x;

z_cg_tfo = A_z;

% Crew: reference point start of cockpit
x_crew = 0.5 * length_cockpit;
x_cg_crew = x_crew + E_x;

z_cg_crew = E_z;

%Fuel: nose of aircraft ****** Important one
s1 = (diam_front / 2)^2 * pi;
s2 = (diam_rear / 2)^2 * pi;
x_fuel = (length_fuselage/4)*(s1 + s2 + 2*(s1*s2)^0.5)/(s1+s2+(s1*s2)^0.5);
x_cg_fuel = x_fuel + A_x;

z_cg_fuel = A_z; % ***********

%Passengers + Baggage: front of cabin: ref front of cabin
x_cabin = 0.5*length_cabin;
x_cg_cabin = x_cabin + F_x;

z_cg_cabin = F_z;

% Assemble into one vector
%x_cg = [x_cg_fuselage x_cg_wing x_cg_canard x_cg_vtail x_cg_vertstab x_cg_engine x_cg_landgear x_cg_fixed_equp x_cg_tfo x_cg_crew x_cg_fuel x_cg_cabin];

x_cg = [x_cg_fuselage x_cg_wing x_cg_canard x_cg_vtail x_cg_engine x_cg_landgear x_cg_fixed_equp x_cg_tfo x_cg_crew x_cg_fuel x_cg_cabin];
z_cg = [z_cg_fuselage z_cg_wing z_cg_canard z_cg_vtail z_cg_engine z_cg_landgear z_cg_fixed_equp z_cg_tfo z_cg_crew z_cg_fuel z_cg_cabin];


%% Determine the CGs for certain scenarios
% Multiply location by weight
temp = x_cg .* Comp_weight;
load1 = sum(temp(1:7))/sum(Comp_weight(1:7)); % Empty weight
load2 = sum(temp(1:9))/sum(Comp_weight(1:9)); % Operating empty weight
load3 = sum(temp(1:10))/sum(Comp_weight(1:10)); % OEW + fuel
load4 = sum(temp(1:11))/sum(Comp_weight(1:11)); % Max take off
load5 = sum(temp([1:9,11]))/sum(Comp_weight([1:9,11]));
CG_motionx = [load1 load2 load3 load4 load5];

temp = z_cg .* Comp_weight;
load1 = sum(temp(1:7))/sum(Comp_weight(1:7)); % Empty weight
load2 = sum(temp(1:9))/sum(Comp_weight(1:9)); % Operating empty weight
load3 = sum(temp(1:10))/sum(Comp_weight(1:10)); % OEW + fuel
load4 = sum(temp(1:11))/sum(Comp_weight(1:11)); % Max take off
load5 = sum(temp([1:9,11]))/sum(Comp_weight([1:9,11]));

CG_motionz = [load1 load2 load3 load4 load5];
weight_motion = [sum(Comp_weight(1:7)) sum(Comp_weight(1:9)) sum(Comp_weight(1:10)) sum(Comp_weight(1:11)) sum(Comp_weight([1:8,10]))];

figure(1) %Plot weight motion and for different load situations
plot(CG_motionx, weight_motion, 'r-', 'Marker', 'o', 'MarkerFaceColor','r')
title('CG Location For loading scenarios')
xlabel('CG Location'); ylabel('Weight')


figure(2) % Plot CG for each load case
scatter(CG_motionx, CG_motionz,'Marker', 'o', 'MarkerFaceColor','k')
title('CG Location in Each Loading Scenario')
xlabel('X Location'); ylabel('Y Location')
axis([0 180 10 30])


figure(3)
scatter(x_cg, z_cg,'Marker', 'o', 'MarkerFaceColor','k')
title('CG Location of each Component')
xlabel('X Location'); ylabel('Y Location')
axis([0 180 10 30])
axis equal

figure(4)
y = zeros(1,length(x_cg));
scatter(x_cg,y, 'Marker', 'o', 'MarkerFaceColor','k')
title('CG Location of Each Component from a Top View')
xlabel('X Location'); ylabel('Y Location')
axis([0 180 10 30])
axis equal