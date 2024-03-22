clc
clear all

r = 5; %mm
l = 200; %mm
T = 1000; %Nm
P = 60000; %N
M = 800; %Nm
I = 0.25*pi*r^4; %mm^4
J = (1/2)*pi*r^4;



sigma_b = (M*r/I); %MPa

sigma_a = P/(pi*r^2); %Mpa

tau_xy = (T*r/J); %Mpa

sigma_x = sigma_b + sigma_a;; %MPa

sigma_y = 0;
sigma_z = 0;
tau_yz = 0;
tau_zx = 0;

I1 = sigma_x + sigma_y + sigma_z;
I2 = (sigma_x*sigma_y + sigma_y*sigma_z + sigma_z*sigma_x - tau_xy^2 - tau_yz^2 - tau_zx^2);
I3 = (sigma_x*sigma_y*sigma_z + 2*tau_xy*tau_yz*tau_zx - sigma_x*tau_yz^2 - sigma_y*tau_yz^2 - sigma_z*tau_xy^2);

sigma123 = roots([1 -I1 I2 -I3])







