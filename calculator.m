%% Fatigue & Fracture Mechanics
%   Term Project
%   Aidan Clark & Patrick Cleary
clear
close all
clc
format shortEng


%% Input
% Prompt user for input
disp('Enter the following parameters:')
% Geometry
radius = input('Shaft Radius (m): ');
length = input('Shaft Length (m): ');
a = input('Snap Ring Depth (m)');

% Material Properties
yield_strength_MPa = input('Material Yield Strength (MPa): ');
k_Ic_MPa    = input('Critical Fracture Toughness Mode I (MPa.m^(1/2)): ');
% k_IIc_MPa   = input('Critical Fracture Toughness Mode II (MPa.m^(1/2)): ');
k_IIIc_MPa  = input('Critical Fracture Toughness Mode III (MPa.m^(1/2)): ');

% Design Constraints
safety_factor_yielding      = input('Safety Factor against Yielding: ');
safety_factor_fracture_I    = input('Safety Factor against Fracture Mode I: ');
% safety_factor_fracture_II   = input('Safety Factor against Fracture Mode II: ');
safety_factor_fracture_III  = input('Safety Factor against Fracture Mode III: ');

% Loading Conditions
torque = input('Applied Torque (N.m): ');
axial_force = input('Axial Force (N): ');
bending_moment = input('Bending Moment (N.m): ');

% Comparison Option
stress_analysis_type = input('Select Comparison Option (1 for Maximum Shear, 2 for Octahedral): ');

% Input Data Unit Conversion
yield_strength = 10e-6 * yield_strength_MPa;
k_Ic = k_Ic_MPa * 10e-6;
% k_IIc = k_IIc_MPa * 10e-6;
k_IIIc = k_IIIc_MPa * 10e-6;



%% Calculations

% Geometric Calculations
alph = a/radius;
beta = 1 - alph;
area = pi * radius^2;
polar_moment_of_inertia = pi/32 * radius^4;

% Torsional stress
tau_xy = torque * radius / polar_moment_of_inertia;

% Axial stress
normal_stress_axial = axial_force / area;

% Bending stress
bending_stress = bending_moment * radius / polar_moment_of_inertia;

% Principal Stress
sigma_1 = 0.5*(normal_stress_axial + bending_stress) + sqrt((0.5*(normal_stress_axial-bending_moment))^2 + tau_xy^2);
sigma_2 = 0.5*(normal_stress_axial + bending_stress) - sqrt((0.5*(normal_stress_axial-bending_moment))^2 + tau_xy^2);

principal_stresses = [sigma_1, sigma_2];

% Effective stress based on user's choice
if stress_analysis_type == 1
    % Maximum shear stress
    effective_stress = max(abs(principal_stresses));
elseif stress_analysis_type == 2
    % Octahedral stress
    effective_stress = sqrt(sum(principal_stresses.^2)/3);
else
    error('Invalid stress analysis type option');
end

% Safety Factor against yielding
safety_factor_yielding_actual = yield_strength / abs(max(principal_stresses));

% Fracture Mechanics
s_g_axial = axial_force/(pi*radius^2);
F_axial = (1/2*beta^1.5)*(1 + 0.5*beta + (3/8)*beta^2 - 0.363*beta^3 + 0.731*beta^4);
k_axial = s_g_axial*F_axial*sqrt(pi*a);

s_g_bending = (4*bending_moment)/(pi*radius^3);
F_bending = (3/(8*beta^2.5))*(1 +0.5*beta + (3/8)*beta^2 + (5/16)*beta^3 + (35/128)*beta^4 + 0.537*beta^5);
k_bending = s_g_bending*F_bending*sqrt(pi*a);

s_g_torsion = (2*torque)/(pi*radius^3);
F_torsion = (3/(8*beta^2.5))*(1 + 0.5*beta + (3/8)*beta^2 + (5/16)*beta^3 + (35/128)*beta^4 + 0.537*beta^5);
k_torsion = s_g_torsion*F_torsion*sqrt(pi*a);


k_I = k_axial + k_bending;
k_III = k_torsion;

% Safety Factors against fracture
safety_factor_fracture_I = k_Ic / k_I;
safety_factor_fracture_III = k_IIIc / k_III;



%% Result Output

% Display results
disp('----- Results -----');
disp('Principal Stresses (Pa):');
disp(principal_stresses);
disp(['Effective Stress (Pa): ', num2str(effective_stress)]);
disp(['Safety Factor against Yielding: ', num2str(safety_factor_yielding_actual)]);
disp(['Safety Factor against Fracture Mode I: ', num2str(safety_factor_fracture_I_actual)]);
disp(['Safety Factor against Fracture Mode II: ', num2str(safety_factor_fracture_II_actual)]);
disp(['Safety Factor against Fracture Mode III: ', num2str(safety_factor_fracture_III_actual)]);

% Write results to a text file
output_file = fopen('output_report.txt', 'w');
fprintf(output_file, '----- Results -----\n');
fprintf(output_file, 'Principal Stresses (Pa):\n');
fprintf(output_file, '%f\n', principal_stresses);
fprintf(output_file, 'Effective Stress (Pa): %f\n', effective_stress);
fprintf(output_file, 'Safety Factor against Yielding: %f\n', safety_factor_yielding_actual);
fprintf(output_file, 'Safety Factor against Fracture Mode I: %f\n', safety_factor_fracture_I);
fprintf(output_file, 'Safety Factor against Fracture Mode III: %f\n', safety_factor_fracture_III);
fclose(output_file);

% Display fracture toughness comparison
disp('----- Fracture Toughness Comparison -----');
disp('Fracture Toughness Values (Pa.m^(1/2):');
disp(fracture_toughness_values);
disp('Critical Stresses for Fracture (Pa):');
disp(critical_stresses);

% Write fracture toughness results to a separate text file
fracture_output_file = fopen('fracture_report.txt', 'w');
fprintf(fracture_output_file, '----- Fracture Toughness Comparison -----\n');
fprintf(fracture_output_file, 'Fracture Toughness Values (Pa.m^(1/2)):\n');
fprintf(fracture_output_file, '%f\n', fracture_toughness_values);
fprintf(fracture_output_file, 'Critical Stresses for Fracture (Pa):\n');
fprintf(fracture_output_file, '%f\n', critical_stresses);
fclose(fracture_output_file);