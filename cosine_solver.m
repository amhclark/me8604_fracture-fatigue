% Aidan Clark, ME8604, amhclark@mun.ca
% Memorial University of Newfoundland
% Faculty of Engineering & Applied Science

% clear;
% close all;
% clc;

global sigma_x sigma_y sigma_z tau_xy tau_yz tau_zx sigma_i
syms l_i m_i n_i

% prompt = {'sigma_x','sigma_y','sigma_z','tau_xy','tau_yz','tau_zx','sigma_i'};
% dlgtitle = 'Input';
% fieldsize = [1 45; 1 45; 1 45; 1 45; 1 45; 1 45; 1 45];
% % definput = {'20','hsv'};
% answer = inputdlg(prompt,dlgtitle,fieldsize);
% 
% sigma_x = str2double(answer{1});
% sigma_y = str2double(answer{2});
% sigma_z = str2double(answer{3});
% tau_xy = str2double(answer{4});
% tau_yz = str2double(answer{5});
% tau_zx = str2double(answer{6});
% sigma_i = str2double(answer{7});

stress = [sigma_x tau_xy tau_zx; tau_xy sigma_y tau_yz; tau_zx tau_yz sigma_z];
[lmn eigvals] = eig(stress);

for k = 1:3
  lmn(:,k) = lmn(:,k)/norm(lmn(:,k));
end

stress*lmn(:,1) - lmn(:,1)*eigvals(1,1);

for i =1:length(eigvals)
    eq1 = (sigma_x - sigma_i)*l_i + tau_xy*m_i + tau_zx*n_i == 0;
    eq2 = tau_xy*l_i + (sigma_y - sigma_i)*m_i + tau_yz*n_i == 0;
    eq3 = tau_zx*l_i + tau_yz*m_i + (sigma_z - sigma_i)*n_i == 0;
    eq4 = l_i^2 + m_i^2 + n_i^2 == 1;

    [A,B] = equationsToMatrix([eq1, eq2, eq3, eq4], [l_i, m_i, n_i]);              
    X = linsolve(A,B);
end

