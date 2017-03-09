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
fuselage = [0.073 0.093 0.122 0.108];
nacelle = [0.017 0.016 0.015];
landing_gear = [0.037 0.051 0.040 0.038];

ave_sizing = [mean(structure) mean(powerplant) mean(fixed_equip) mean(empty_weight) mean(wing) mean(empennage) mean(fuselage) mean(nacelle) mean(landing_gear)];
weight_val = ave_sizing .* GW;

Weight_frac = [mean(fuselage) mean(wing) mean(empennage) mean(powerplant) mean(landing_gear) mean(fixed_equip) 0.0028 0.0033 0.456 0.036];
Comp_weight = Weight_frac .* GW;

%% Find the overall locations of each component
    % Start with reference locations
A_x = 20; % Nose
B_x = 120; % front tip of wing
C_x = 180; % front of vert stab
D_x = 170; % engine inlet
E_x = 30; % front of cockpit
F_x = 45; % front of cabin

%% Componenet cg locations
    % all starred quantities must be fillled in from geometry
    
% Fueslage: reference point, nose of aircraft
length_fuselage = 150;%**** in feet for twin engine aft mounted fuselage
x_fuselage = 0.485 * length_fuselage;
x_cg_fuselage = x_fuselage + A_x;

% wing: reference point is the leading tip of the swept wing
c_root = 50;%*** Root chord
le_sweep = 40; % *** Leading edge sweep angle
te_sweep = 20; % *** Trailing edge sweep
b = 100;% **** Wing span b
x_wing = 0.70*(1-0.35)*(b/2)*(tand(le_sweep)-tand(te_sweep)) + (tand(le_sweep)*(0.35*(b/2)));
x_cg_wing = x_wing + B_x;

% Empennage: reference point: leading edge of tail root 
c_root_tail = 10;%*** Root chord of vertical stabilizer
vle_sweep = 40; % *** Leading edge sweep angle of vertical stablizer
vte_sweep = 20; % *** Trailing edge sweep of vertical stabilizer
v_b = 15;% **** vertical tail length
x_vertstab = 0.42*(1-0.38)*v_b*(tand(vle_sweep)-tand(vte_sweep)) + (tand(vle_sweep)*(0.38*(v_b))); % for no horizontal sabilizer
%x_cg_vert_stab = 0.42*(1-0.55)*v_b*(tand(vle_sweep)-tand(vte_sweep)) +
    %(tand(vle_sweep)*(0.55*(v_b))) % if T tail configuration
x_cg_vertstab = x_vertstab + C_x;
    
%Engine group : reference point inlet of engine **** needs revision
length_engine = 14.05; % Length of engine from P&W specs
x_engine = 0.5*length_engine;
x_cg_engine = x_engine + D_x;

%Landing gear: reference point, nose of the aircraft
    %assume mostly vertical struts so cg is where the landing gear are
    %placed
x_landgear = 76; %****
x_cg_landgear = x_landgear + A_x;

% Fixed Equipment: ref point nose of aircraft
x_fixed_equp = 0.3 * length_fuselage; % *** Needs Far Mosre Detail
x_cg_fixed_equp = x_fixed_equp + A_x;

% Trapped fuel and oil: ref point, nose of aircraft
    % near engine case
x_tfo = 0; %**** assuming its near engine inlet
x_cg_tfo = x_tfo + D_x;

% Crew: reference point start of cockpit
length_cockpit = 10; %****
x_crew = 0.5 * length_cockpit;
x_cg_crew = x_crew + E_x;

%Fuel: nose of aircraft ****** Important one
diam_front = 13; % **** diameter of fuselage at front 
diam_rear = 9; % Diameter of fuselage at rear
s1 = (diam_front / 2)^2 * pi;
s2 = (diam_rear / 2)^2 * pi;
x_fuel = (length_fuselage/4)*(s1 + s2 + 2*(s1*s2)^0.5)/(s1+s2+(s1*s2)^0.5);
x_cg_fuel = x_fuel + A_x;

%Passengers + Baggage: front of cabin: ref front of cabin
length_cabin = 75;% *****
x_cabin = 0.5*length_cabin;
x_cg_cabin = x_cabin + F_x;

% Assemble into one vector
x_cg = [x_cg_fuselage x_cg_wing x_cg_vertstab x_cg_engine x_cg_landgear x_cg_fixed_equp x_cg_tfo x_cg_crew x_cg_fuel x_cg_cabin];

%% Determine the CGs for certain scenarios
% Multiply location by weight
temp = x_cg .* Comp_weight;
load1 = sum(temp(1:6))/sum(Comp_weight(1:6)) % Empty weight
load2 = sum(temp(1:8))/sum(Comp_weight(1:8)) % Operating empty weight
load3 = sum(temp(1:9))/sum(Comp_weight(1:9)) % OEW + fuel
load4 = sum(temp(1:10))/sum(Comp_weight(1:10)) % Max take off
load5 = sum(temp([1:8,10]))/sum(Comp_weight([1:8,10]))

CG_motion = [load1 load2 load3 load4 load5];
weight_motion = [sum(Comp_weight(1:6)) sum(Comp_weight(1:8)) sum(Comp_weight(1:9)) sum(Comp_weight(1:10)) sum(Comp_weight([1:8,10]))]

plot(CG_motion, weight_motion, 'r-', 'Marker', 'o', 'MarkerFaceColor','g')
title('CG Location For loading scenarios')
xlabel('CG Location'); ylabel('Weight')

