%% Fatigue & Fracture Mechanics
%   Term Project
%   Aidan Clark & Patrick Cleary
%   ME 8604

clear
close all
clc
format shortEng

%% Input
% Prompt user for input
disp('Enter the following parameters:')
% Geometry Input
radius  = input('Shaft Radius (m): ');
length  = input('Shaft Length (m): ');
a       = input('Snap Ring Depth (m): ');

% Material Properties Input
yield_strength_MPa  = input('Material Yield Strength (MPa): ');
k_Ic                = input('Critical Fracture Toughness Mode I (MPa.m^(1/2)): ');
k_IIIc              = input('Critical Fracture Toughness Mode III (MPa.m^(1/2)): ');

% Design Constraints
safety_factor_yielding_design      = input('Safety Factor against Yielding: ');
safety_factor_fracture_I_design    = input('Safety Factor against Fracture Mode I: ');
safety_factor_fracture_III_design  = input('Safety Factor against Fracture Mode III: ');

% Loading Conditions
torque          = input('Applied Torque (N.m): ');
axial_force     = input('Axial Force (N): ');
bending_moment  = input('Bending Moment (N.m): ');

% Comparison Option
stress_analysis_type = input('Select Comparison Option (1 for Maximum Shear, 2 for Octahedral): ');

% Input Data Unit Conversion (for future calculations)
yield_strength = 1e6 * yield_strength_MPa; 


%% Part I Calculations
% Geometric Calculations
area = pi * radius^2;                       % in m^2
polar_moment_of_inertia = (pi*radius^4)/2;  % in m^4
moment_of_inertia = (pi*radius^4)/4;        % in m^4

% Torsional stress
torsional_stress = torque * radius / polar_moment_of_inertia;

% Axial stress
normal_stress_axial = axial_force / area;

% Bending stress
bending_stress = bending_moment * radius / moment_of_inertia;

% Principal Stress Setup
sigma_x = bending_stress + normal_stress_axial;
sigma_y = 0;
sigma_z = 0;
tau_xy  = torsional_stress;
tau_yz  = 0;
tau_zx  = 0;

% 3D stress state solution
I1 = sigma_x + sigma_y + sigma_z;
I2 = sigma_x*sigma_y + sigma_y*sigma_z + sigma_z*sigma_x - tau_xy^2 - tau_yz^2 - tau_zx^2;
I3 = sigma_x*sigma_y*sigma_z + 2*tau_xy*tau_yz*tau_zx - sigma_x*tau_yz^2 - sigma_y*tau_zx^2 - sigma_z*tau_xy^2;

% Solving Cubic principal stress equation
principal_stresses  = sort(roots([1 -1*I1 +I2 -I3]));    % Sorted roots of cubic 3d stress equation
sigma_1_MPa         = (principal_stresses(3))*1e-6;            % Principal stress definition based on common practice nomenclature
sigma_2_MPa         = (principal_stresses(1))*1e-6;
sigma_3_MPa         = (principal_stresses(2))*1e-6;

% New list defined for principal stresses in MPa
principal_stresses_MPa = [sigma_1_MPa, sigma_2_MPa, sigma_3_MPa];

% Effective stress based on user's choice
if stress_analysis_type == 1
    % Maximum shear stress
    effective_stress = max([abs(sigma_1_MPa - sigma_2_MPa) abs(sigma_2_MPa - sigma_3_MPa) abs(sigma_3_MPa - sigma_1_MPa)]);
elseif stress_analysis_type == 2
    % Octahedral stress
    effective_stress = (1/sqrt(2))*sqrt((sigma_1_MPa - sigma_2_MPa)^2 + (sigma_2_MPa - sigma_3_MPa)^2 + (sigma_3_MPa - sigma_1_MPa)^2);
else
    error('Invalid stress analysis type option');
end

% Safety Factor
safety_factor = yield_strength_MPa/effective_stress;


%% Part II Calculations
% Geometric Features
Alpha = a/radius; Beta = 1 - Alpha; h = length/2;

% Fracture Mechanics (S_g F and K for axial, bending and torsion)
s_g_axial = (axial_force/(pi*radius^2))*1e-6;  % MPa
F_axial = (1/(2*Beta^(1.5)))*(1 + 0.5*Beta + (3/8)*Beta^2 - 0.363*Beta^3 + 0.731*Beta^4);
k_axial = s_g_axial*F_axial*sqrt(pi*a);         % MPa sqrt(m)

s_g_bending = ((4*bending_moment)/(pi*radius^3))*1e-6;     % MPa
F_bending = (3/(8*Beta^(2.5)))*(1 + 0.5*Beta + (3/8)*Beta^2 + (5/16)*Beta^3 + (35/128)*Beta^4 + 0.537*Beta^5);
k_bending = s_g_bending*F_bending*sqrt(pi*a);               % MPa sqrt(m)

s_g_torsion = ((2*torque)/(pi*radius^3))*1e-6; % MPa
F_torsion = (3/(8*Beta^2.5))*(1 + 0.5*Beta + (3/8)*Beta^2 + (5/16)*Beta^3 + (35/128)*Beta^4 + 0.208*Beta^5);
k_torsion = s_g_torsion*F_torsion*sqrt(pi*a);   % MPa sqrt(m)

% Calculation of fracture toughness in Mode I and Mode III
k_I = k_axial + k_bending;
k_III = k_torsion;

% LEFM Check
LEFM_kI = (4/pi)*(k_I/yield_strength_MPa)^2;
LEFM_kIII = (4/pi)*(k_III/yield_strength_MPa)^2;

if a >= LEFM_kI && a >= LEFM_kIII && (radius - a) >= LEFM_kI && (radius - a) >= LEFM_kIII && h >= LEFM_kI && h >= LEFM_kIII
    LEFM_valid = 'True';
else 
    LEFM_valid = 'False';
end

% Plane Strain Check
PS_kI = 2.5*(k_I/yield_strength_MPa)^2;
PS_kIII = 2.5*(k_III/yield_strength_MPa)^2;

if a >= PS_kI && a >= PS_kIII && (radius - a) >= PS_kI && (radius - a) >= PS_kIII
    PS_valid = 'True';
else
    PS_valid = 'False';
end

% Safety Factors against fracture
safety_factor_fracture_axial = k_Ic / k_axial;
safety_factor_fracture_bending = k_Ic / k_bending;
safety_factor_fracture_torsion = k_IIIc / k_torsion;
safety_factor_fracture_I = k_Ic / k_I;
safety_factor_fracture_III = k_IIIc / k_III;

%Torsional Yield wrt Octahedral Stress relationship
torsion_yield = yield_strength_MPa/sqrt(3); 
torsion_yield_Pa = torsion_yield*1e6;

plastic_force = yield_strength*pi*radius^2*(1-Alpha)^2;
plastic_moment = (4/3)*radius^3*yield_strength*(1-Alpha)^3; 
plastic_torque = (2/3)*pi*radius^3*(1-Alpha)^3*torsion_yield_Pa;

% Safety Factors against fully plastic yielding
plastic_force_fos = plastic_force/axial_force;
plastic_moment_fos = plastic_moment/bending_moment;
plastic_torque_fos = plastic_torque/torque;


%% FOS Comparisons (actual vs design)
% Comparing User Input FOS vs actual FOS to determine if the loading
% conditions allow for safe operation

% Yielding Comparison
if safety_factor >= safety_factor_yielding_design
    yield_comparison = 'Safe';
else
    yield_comparison = 'Unsafe';
end

% Mode I Fracture Comparison
if safety_factor_fracture_I >= safety_factor_fracture_I_design
    FOS_modeI_comparison = 'Safe';
else
    FOS_modeI_comparison = 'Unsafe';
end

% Mode III Fracture Comparison
if safety_factor_fracture_III >= safety_factor_fracture_III_design
    FOS_modeIII_comparison = 'Safe';
else
    FOS_modeIII_comparison = 'Unsafe';
end


%% Result Output

% Specifying precision for .txt file output
precision = 3;
formatSpec = ['%.', num2str(precision), 'f'];

% Part I: Results .txt File
output_file_part_I = fopen('un-notched_mechanics_output.txt', 'w');
fprintf(output_file_part_I,'----- User Input -----\n');
fprintf(output_file_part_I, ['Radius (m) ', formatSpec, '\n'], radius);
fprintf(output_file_part_I, ['Length (m) ', formatSpec, '\n'], length);
fprintf(output_file_part_I, ['Yield Strength (MPa) ', formatSpec, '\n'], yield_strength_MPa);
fprintf(output_file_part_I, ['Torque (Nm) ', formatSpec, '\n'], torque);
fprintf(output_file_part_I, ['Axial Force (N) ', formatSpec, '\n'], axial_force);
fprintf(output_file_part_I, ['Bending Moment (N/m) ', formatSpec, '\n'], bending_moment);
fprintf(output_file_part_I, ' \n');
fprintf(output_file_part_I, '----- Results -----\n');
fprintf(output_file_part_I, 'Principal Stresses (MPa):\n');
fprintf(output_file_part_I, [formatSpec, '\n'], principal_stresses_MPa);
fprintf(output_file_part_I, 'Stress Analysis Type (1 = Max Shear, 2 = Octahedral): %s\n',num2str(stress_analysis_type));
fprintf(output_file_part_I, ['Effective Stress (MPa): ', formatSpec, '\n'], effective_stress);
fprintf(output_file_part_I, ['Safety Factor: ', formatSpec, '\n'], safety_factor);
fprintf(output_file_part_I, 'Is the yield FOS safe/unsafe compared with design values: %s\n', yield_comparison);
fclose(output_file_part_I);
disp('Report of Un-Notched Shaft Analysis (Part I) saved as: un-notched_mechanics_output.txt')

% Part II: Results .txt File
output_file_part_II = fopen('fracture_mechanics_output.txt', 'w');
fprintf(output_file_part_II,'----- User Input -----\n');
fprintf(output_file_part_II, ['Radius (m) ', formatSpec, '\n'], radius);
fprintf(output_file_part_II, ['Length (m) ', formatSpec, '\n'], length);
fprintf(output_file_part_II, ['Yield Strength (MPa) ', formatSpec, '\n'], yield_strength_MPa);
fprintf(output_file_part_II, ['Torque (Nm) ', formatSpec, '\n'], torque);
fprintf(output_file_part_II, ['Axial Force (N) ', formatSpec, '\n'], axial_force);
fprintf(output_file_part_II, ['Bending Moment (N/m) ', formatSpec, '\n'], bending_moment);
fprintf(output_file_part_II, ['Snap Ring Depth (m): ', formatSpec, '\n'], a);
fprintf(output_file_part_II, ['Fracture Toughness I (MPa', char(8730),'m): ', formatSpec, '\n'], k_Ic);
fprintf(output_file_part_II, ['Fracture Toughness III (MPa', char(8730),'m): ', formatSpec, '\n'], k_IIIc);
fprintf(output_file_part_II, ['Safety Factor Yielding (Design): ', formatSpec, '\n'], safety_factor_yielding_design);
fprintf(output_file_part_II, ['Safety Factor Fracture I (Design): ', formatSpec, '\n'], safety_factor_fracture_I_design);
fprintf(output_file_part_II, ['Safety Factor Fracture III (Design): ', formatSpec, '\n'], safety_factor_fracture_III_design);
fprintf(output_file_part_II, ' \n');
fprintf(output_file_part_II, '----- Results -----\n');
fprintf(output_file_part_II, 'LEFM Valid: %s\n', LEFM_valid);
fprintf(output_file_part_II, 'Plane Strain Valid: %s\n', PS_valid);
fprintf(output_file_part_II, ['Fracture Toughness, Axial (MPa', char(8730),'m): ', formatSpec, '\n'], k_axial);
fprintf(output_file_part_II, ['Fracture Toughness, Bending (MPa', char(8730),'m): ', formatSpec, '\n'], k_bending);
fprintf(output_file_part_II, ['Fracture Toughness, Torsion (MPa', char(8730),'m): ', formatSpec, '\n'], k_torsion);
fprintf(output_file_part_II, ['Fracture Toughness, Mode I (MPa', char(8730),'m): ', formatSpec, '\n'], k_I);
fprintf(output_file_part_II, ['Fracture Toughness, Mode III (MPa', char(8730),'m): ', formatSpec, '\n'], k_III);
fprintf(output_file_part_II, ['Torsion Yield Strength (MPa): ', formatSpec, '\n'], torsion_yield);
fprintf(output_file_part_II, ['Plastic Force, Axial (N): ', formatSpec, '\n'], plastic_force);
fprintf(output_file_part_II, ['Plastic Force, Bending (Nm): ', formatSpec, '\n'], plastic_moment);
fprintf(output_file_part_II, ['Plastic Force, Torsion (Nm): ', formatSpec, '\n'], plastic_torque);
fprintf(output_file_part_II, ['Fracture FOS (Axial): ', formatSpec, '\n'], safety_factor_fracture_axial);
fprintf(output_file_part_II, ['Fracture FOS (Bending): ', formatSpec, '\n'], safety_factor_fracture_bending);
fprintf(output_file_part_II, ['Fracture FOS (Torsion): ', formatSpec, '\n'], safety_factor_fracture_torsion);
fprintf(output_file_part_II, ['Fracture FOS (Mode I): ', formatSpec, '\n'], safety_factor_fracture_I);
fprintf(output_file_part_II, ['Fracture FOS (Mode III): ', formatSpec, '\n'], safety_factor_fracture_III);
fprintf(output_file_part_II, ['Plastic Yield FOS (Axial): ', formatSpec, '\n'], plastic_force_fos);
fprintf(output_file_part_II, ['Plastic Yield FOS (Bending): ', formatSpec, '\n'], plastic_moment_fos);
fprintf(output_file_part_II, ['Plastic Yield FOS (Torsion): ', formatSpec, '\n'], plastic_torque_fos);
fprintf(output_file_part_II, 'Is the Mode I Fracture FOS safe/unsafe compared with design value: %s\n', sprintf(FOS_modeI_comparison));
fprintf(output_file_part_II, 'Is the Mode III Fracture FOS safe/unsafe compared with design value: %s\n', FOS_modeIII_comparison);

fclose(output_file_part_II);
disp('Report of Fracture Analysis (Part II) saved as: fracture_mechanics_output.txt')
