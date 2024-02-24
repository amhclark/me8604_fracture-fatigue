% User Input
disp('Please enter the following information:');

% Geometry of the shaft
radius = input('Shaft Radius (m): ');
length = input('Shaft Length (m): ');

% Material properties
yield_strength = input('Yield Strength (Pa): ');
fracture_toughness_mode_I = input('Critical Fracture Toughness (Mode I) (Pa.m^0.5): ');
fracture_toughness_mode_II = input('Critical Fracture Toughness (Mode II) (Pa.m^0.5): ');
fracture_toughness_mode_III = input('Critical Fracture Toughness (Mode III) (Pa.m^0.5): ');

% Design constraints
safety_factor_yielding = input('Safety Factor Against Yielding: ');
safety_factor_mode_I = input('Safety Factor Against Fracture (Mode I): ');
safety_factor_mode_II = input('Safety Factor Against Fracture (Mode II): ');
safety_factor_mode_III = input('Safety Factor Against Fracture (Mode III): ');

% Loading conditions
torque = input('Torque (N.m): ');
axial_force = input('Axial Force (N): ');
bending_moment = input('Bending Moment (N.m): ');

% Loading option
loading_option = input('Choose loading option (1 for Maximum Shear, 2 for Octahedral): ');

% Calculate three-dimensional stresses
% (You may replace the following with your own stress calculation function)
% Assuming simple solid shaft with circular cross-section
shear_stress = (torque * radius) / (pi * radius^3 / 2);
axial_stress = axial_force / (pi * radius^2);
bending_stress = bending_moment / (pi * radius^3 / 2);

% Principal stresses
principal_stresses = eig([shear_stress, 0, 0; 0, axial_stress, 0; 0, 0, bending_stress]);

% Evaluate effective stress based on user input
if loading_option == 1
    effective_stress = sqrt(principal_stresses(1)^2 + principal_stresses(2)^2 + principal_stresses(3)^2);
elseif loading_option == 2
    effective_stress = sqrt((principal_stresses(1) - principal_stresses(2))^2 + (principal_stresses(2) - principal_stresses(3))^2 + (principal_stresses(3) - principal_stresses(1))^2);
else
    error('Invalid loading option. Please choose 1 or 2.');
end

% Evaluate safety factor
safety_factor = min([yield_strength / safety_factor_yielding, ...
                     fracture_toughness_mode_I / safety_factor_mode_I, ...
                     fracture_toughness_mode_II / safety_factor_mode_II, ...
                     fracture_toughness_mode_III / safety_factor_mode_III]);

% Print the output in a text file
output_file = fopen('shaft_analysis_report.txt', 'w');
fprintf(output_file, 'Shaft Analysis Report:\n\n');
fprintf(output_file, 'Geometry:\n');
fprintf(output_file, '  Radius: %.4f m\n', radius);
fprintf(output_file, '  Length: %.4f m\n\n', length);

fprintf(output_file, 'Material Properties:\n');
fprintf(output_file, '  Yield Strength: %.4f Pa\n', yield_strength);
fprintf(output_file, '  Critical Fracture Toughness (Mode I): %.4f Pa.m^0.5\n', fracture_toughness_mode_I);
fprintf(output_file, '  Critical Fracture Toughness (Mode II): %.4f Pa.m^0.5\n', fracture_toughness_mode_II);
fprintf(output_file, '  Critical Fracture Toughness (Mode III): %.4f Pa.m^0.5\n\n', fracture_toughness_mode_III);

fprintf(output_file, 'Design Constraints:\n');
fprintf(output_file, '  Safety Factor Against Yielding: %.4f\n', safety_factor_yielding);
fprintf(output_file, '  Safety Factor Against Fracture (Mode I): %.4f\n', safety_factor_mode_I);
fprintf(output_file, '  Safety Factor Against Fracture (Mode II): %.4f\n', safety_factor_mode_II);
fprintf(output_file, '  Safety Factor Against Fracture (Mode III): %.4f\n\n', safety_factor_mode_III);

fprintf(output_file, 'Loading Conditions:\n');
fprintf(output_file, '  Torque: %.4f N.m\n', torque);
fprintf(output_file, '  Axial Force: %.4f N\n', axial_force);
fprintf(output_file, '  Bending Moment: %.4f N.m\n\n', bending_moment);

fprintf(output_file, 'Loading Option: %d (1 for Maximum Shear, 2 for Octahedral)\n\n', loading_option);

fprintf(output_file, 'Results:\n');
fprintf(output_file, '  Three-Dimensional Stresses:\n');
fprintf(output_file, '    Shear Stress: %.4f Pa\n', shear_stress);
fprintf(output_file, '    Axial Stress: %.4f Pa\n', axial_stress);
fprintf(output_file, '    Bending Stress: %.4f Pa\n\n', bending_stress);

fprintf(output_file, '  Principal Stresses:\n');
fprintf(output_file, '    Sigma1: %.4f Pa\n', principal_stresses(1));
fprintf(output_file, '    Sigma2: %.4f Pa\n', principal_stresses(2));
fprintf(output_file, '    Sigma3: %.4f Pa\n\n', principal_stresses(3));

fprintf(output_file, '  Effective Stress: %.4f Pa\n\n', effective_stress);

fprintf(output_file, '  Safety Factor: %.4f\n', safety_factor);

fclose(output_file);

disp('Analysis completed. Results saved in "shaft_analysis_report.txt".');
