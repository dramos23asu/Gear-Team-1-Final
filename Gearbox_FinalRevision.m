clc; clear ; close all 

% Created By: Eric Berkely, David Ramos, Jose Segovia, Cameron Shankland
% Contact through Slack for questions/concerns.

%% Instructions

% Before running this code, make sure that you have downloaded the folder
% together. Do not just run the m-file by itself as there are some files
% from the folder needed to run this code. 

% 1) Input meterial properties

% 2) Input range of gear box heights you wish to test

% 3) Run the code

% 4) Output plots show safety factors vs gear box height of critical gears,
% in this case pinion 3 since it is the smallest gear transmitting the
% largest load

% 5) Find output table with closest critical safety factor, this table gives all
% important information

% 6) If there is no output table, check input values or material
% properties as there may be no valid solution. You can also check the
% plots to see where the safety factors are at for your inputs. 

%% Input Parameters (Variable)

% Meterial Properties (current: Steel, Flame of Induction Hardened) 

% Allowable contact stress number 
Sc2 = 170000; % Gears 2 and 4 
Sc3 = 170000; % Pinions 3 and 5

% Allowable bending stress number (psi) (table 14-3)
St_5 = 50000; % Pinion 5
St_4 = 50000; % Gear 4 
St_3 = 50000; % Pinion 3 
St_2 = 50000; % Gear 2

% Elastic Coefficient 
Cp = 2300;

% Gearbox Dimensions (input range of heights)
gearbox_height_vec = 5:.1:10; 

% Keep these constant
clearance = 1; % in
case_thickness = 1; % in

%% Input Parameters (Constant)

% Allowable gear train values
train_min = 17;
train_max = 21;
gear_train = (train_min:1:train_max);
gr = 18; %nominal gear train value

% Input Torque
T_in = 38.6; % ft-lb

% Input Speed
w_in = 200; % RPM

% Safety Factors
surface_n = 1.2; 
bending_n = 1.5;

%% Optimization
for i = 1:length(gearbox_height_vec)

gearbox_height(i) = gearbox_height_vec(1,i);

% Pressure Angle
p_angle = 20; % degrees
p_rad = 20*(pi/180); % converting to radians
k = 1; % For full depth teeth
mg = sqrt(gr); % Ratio of mating gears

% Minimum Number of Teeth
N_min = ceil((2*k)/((1+2*mg)*(sin(p_rad))^2)*(mg+sqrt(mg^2+(1+2*mg)*(sin(p_rad))^2))); % Eq. 13-11
        % Always rounding up

% Number of Teeth for Each Gear
        % Assuming that each alternate gear is equal to each other

% Number of Teeth for Pinions 
N3 = N_min;
N5 = N_min;

% Number of Teeth for Gears
N2 = floor(sqrt(gr)*N3);
N4 = N2;
        % Always rounding down (to maintain gear ratio)

P_min_exact = (N3+(N2/2)+(N5/2)+2)/(gearbox_height(i)- ...
    (clearance+case_thickness));
P_min(i) = ceil(P_min_exact);
        % Always rounding up

% Face Width (3-5) Times Circular Pitch
F(i) = round(4*(pi/P_min(i)),1);

% Resulting Speeds for Each Gear/Pinion (RPM)
w2 = w_in;
w3 = (N2/N3)*w_in;
w4 = w3;
w5 = (N2/N3)^2*w_in;

% Diametral Pitches for each Gear/Pinion
d2(i) = N2/P_min(i);
d4(i) = d2(i);
d3(i) = N3/P_min(i);
d5(i) = d3(i);

% Pitch-Line Velocities
V23 = (pi*d2(i)*w2)/12; % ft/min
V45 = (pi*d5(i)*w5)/12; % ft/min

% Power
H = T_in*w_in/5252; % hp

% Torques on each Gear/Pinion
T2 = T_in; %ft-lb
T3 = T2*(w2/w3); % ft-lb
T4 = T3; % ft-lb
T5 = T2*(w2/w5); % ft-lb

% Transmitted Loads (Tangential)
W23 = 33000*(H/V23); % loads between 2/3
W45 = 33000*(H/V45); % loads between 4/5

%% Gear 4 Wear
I = ((cos(p_rad)*sin(p_rad))/2)*(sqrt(gr)/(sqrt(gr)-1));

Qv = 10; % Quality Number (Precision)
B = 0.25*(12-Qv)^(2/3);
A = 50+56*(1-B);
Kv = ((A+sqrt(V45))/A)^(B);

% Load Distribution Factor Constants

if F(i) <= 1
    Cpf = (F(i)/(10*d3(i)))-0.025;
else
    Cpf = (F(i)/(10*d3(i)))-0.0375+0.0125*F(i);
end

Cmc = 1; % uncrowned teeth

Cpm = 1; % may depend on location of where the gears are placed on the shaft, assume 1

% Table 14-9, Precision Enclosed Units

Cma = 0.0675 + 0.0128*F(i) + (-0.926*10^-4)*(F(i)^2);

Ce = 1; % all other conditions

% Load Distribution Factor, Km
Km = 1 + Cmc*(Cpf*Cpm+Cma*Ce);

% AGMA Stress 
sigma_c_4 = Cp*sqrt((W45*Kv*Km)/(d5(i)*F(i)*I));

% Life Factor
operation_time = 10; % years
L4 = (operation_time*8766)*60*w4; % total cycles
Zn = 1.4488*L4^-0.023; % Figure 14-15

Kr = 1; % Reliability Factor
Kt = 1; % Temperature Factor
Ch = 1; % Hardness ratio factor

Sc_4 = (surface_n*sigma_c_4)/Zn;

n_c4(i) = (Sc2*Zn)/sigma_c_4; %gear 4 wear safety factor 

%% Gear 4 Bending

J_4 = 0.39; % Geometry Factor
KB_4 = 1;

sigma_b_4 = (W45*Kv*P_min(i)*Km)/(F(i)*J_4); % psi

Yn_4 = 1.3558*L4^-0.0178; % Figure 14-14

sigma_all_4 = Yn_4*St_4;

n_b4(i) = sigma_all_4/sigma_b_4; % gear for bending safety factor

%% Gear 5 Bending and Wear

J_5 = 0.26; % changes with gear teeth

L_5 = (operation_time*8766)*60*w5;

Yn_5 = 1.3558*L_5^-0.0178; % Figure 14-14

Zn_5 = 1.4488*L_5^-0.023; % Figure 14-15

% Contact Stress on 5
sigma_c_5 = Cp*sqrt((W45*Kv*Km)/(d5(i)*F(i)*I));

% Bending Stress on 5
sigma_b_5 = (W45*Kv*P_min(i)*Km*KB_4)/(F(i)*J_5);

Sc_5 = (surface_n*sigma_c_5)/Zn_5;

sigma_all_5 = Yn_5*St_5; 

% Resulting Safety Factors
n_c5(i) = (Sc3*Zn_5)/sigma_c_5; % wear safety factor
n_b5(i) = sigma_all_5/sigma_b_5; % bending safety factor

%% Gear 2 Bending and Wear 

J_2 = J_4;

kv_2 = ((A+sqrt(V23))/A)^(B);
sigma_c_2 = Cp*sqrt((W23*kv_2*Km)/(d3(i)*F(i)*I));

L2 = (operation_time*8766)*60*w2;
Zn_2 = 1.4488*L2^-0.023;
Sc_2 = (surface_n*sigma_c_2)/Zn_2;
n_c2(i) = (Sc2*Zn_2)/sigma_c_2; % gear two safety factor wear

Yn_2 = 1.3558*L2^-0.0178; % Figure 14-14

sigma_b_2 = (W23*kv_2*P_min(i)*Km)/(F(i)*J_2);
sigma_all_2 = St_2*Yn_2;
n_b2(i) = sigma_all_2/sigma_b_2; %gear 2 bending safety factor

%% Gear 3 Bending and Wear

L_3 = (operation_time*8766)*60*w3;
J_3 = J_5;
Yn_3 = Yn_4;
Zn_3 = Zn; 

sigma_c_3 = Cp*sqrt((W23*kv_2*Km)/(d3(i)*F(i)*I)); % wear stress

sigma_b_3 = (W23*kv_2*P_min(i)*Km*KB_4)/(F(i)*J_3); % bending stress

sigma_all_3 = Yn_3*St_3;

Sc_3 = (surface_n*sigma_c_3)/Zn_3;

% Safety Factors
n_c3(i) = (Sc3*Zn_3)/sigma_c_3; % wear 
n_b3(i) = sigma_all_3/sigma_b_3; % bending


%% Output Tables 

if n_c3(i) >= surface_n && n_b3(i) >= bending_n

optimal_location(i) = i;

else

end

end

opt_index = min(optimal_location (optimal_location > 0));

i = opt_index;
Gears = ["Gear 2";"Gear 4";"Pinion 3";"Pinion 5"];
Diameters_in = [d2(i);d4(i);d3(i);d5(i)];
Teeth_Number = [N2;N4;N3;N5];
Face_Width_in = [F(i);F(i);F(i);F(i)];
SF_Bending = [n_b2(i);n_b4(i);n_b3(i);n_b5(i)];
SF_Wear = [n_c2(i);n_c4(i);n_c3(i);n_c5(i)];
Rotational_Speed_RPM = [w2;w4;w2;w5];
Tangential_Load_lbf = [W23;W45;W23;W45];
Diametral_Pitch = [P_min(i);P_min(i);P_min(i);P_min(i)];
Gearbox_Height_in = gearbox_height(i);
% Display best optimized design based on safety factors

Optimized_Output = table(Gears,Diameters_in,Teeth_Number,Diametral_Pitch,Face_Width_in,SF_Wear, ...
    SF_Bending,Rotational_Speed_RPM,Tangential_Load_lbf)

fprintf('Optimized Gearbox Height = %2.1f in\n',Gearbox_Height_in)

%% Safety Factor vs Gearbox Height Plots

% Drawing horizontal line for safety factors
x_b = [gearbox_height_vec(1) gearbox_height_vec(length(gearbox_height_vec))];
y_b = [bending_n bending_n];
y_c = [surface_n surface_n];

% Image of Resulting Design
gb_design = imread('WindMillGearboxDesign.jpg');
figure
imshow(gb_design)

figure
hold on
plot(gearbox_height,n_c3)
plot(gearbox_height,n_b3)
line(x_b,y_b, 'Color', 'r') % Bending safety factor line
line(x_b, y_c, 'Color', 'b') % Surface safety factor line
title('Gear 3 Safety Factors Vs Gearbox Height')
xlabel('Gear Box Height (in)')
ylabel('Safety Factor Value')

legend('Surface Safety Factor','Bending Safety Factor','Minimum Bending Safety Factor','Minimum Surface Safety Factor')

