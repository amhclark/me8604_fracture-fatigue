function result = stress(radius, axial_force, bending_moment, torque)
    cross_section = 2*pi*radius^2;
    polar_moment = (pi*radius^4)/2;
    i_z = (pi*radius^4)/4;

    axial_stress = axial_force/cross_section;
    torque_shear = (torque*radius)/polar_moment;
    bending_stress = (bending_moment)/i_z;
    principal_stress = sqrt(axial_stress^2 + bending_stress^2 + 3*torque_shear^2);

    result = principal_stress;
end
