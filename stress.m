function result = stress(radius, axial_force, bending_moment, torque)
    
    cross_section = 2*pi*radius^2;
    polar_moment = (pi*radius^4)/2;
    i_z = (pi*radius^4)/4;

    axial_stress = axial_force/cross_section;
    torque_shear = (torque*radius)/polar_moment;
    bending_stress = (bending_moment)/i_z;

    % TODO: replace with real principal calculation (placeholder for debug)
    principal_stress = axial_stress+torque_shear+bending_stress;
    result = principal_stress;

end


function result = effective_stress(principal_stress, stress_type)

    % TODO: calculation step to isolate, sig1 sig2 sig3 for calclation

    if stress_type == vonmises
        eff_stress = (1/sqrt(2))*sqrt((sigma_1 - sigma_2)^2 + (sigma_2 - sigma_3)^2 + (sigma_3 - sigma_1)^2);
        result = eff_stress;
    end

    % TODO: replace if with elif
    if stress_type == tresca
        eff_stress = max(abs(sigma_1 - sigma_2), abs(sigma_2 - sigma_3), abs(sigma_3 - sigma_1));
        result = eff_stress;
    end    

end

% TODO: write function
function result = safety(effective_stress, safety_yield, safety_fracture)

end
