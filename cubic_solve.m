clear;
close all;
clc;

% Aidan Clark, ME8604, amhclark@mun.ca
% Memorial University of Newfoundland
% Faculty of Engineering & Applied Science

F = [1 -90 -4844 189240];       % Coefficients of cubic equation from big to small
roots(F)                        % Solve the roots of the equation


% Solve Graphically
j = 1;
for i = -150:150
    
    RT(j) = i^3 - 90*i^2 - 4844*i +189240;
    j = j+1;
    xt = -150:150;

end

plot(xt,RT)
xline(0);
yline(0);

