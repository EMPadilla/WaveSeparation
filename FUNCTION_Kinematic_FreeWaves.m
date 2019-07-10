%% FUNCTION_Kinematic_FreeWaves
%  This function provides the theoretical Kinematics for free waves: 
%         -- Wavenumber (k)
%         -- Celerity (c), 
%         -- Group Celerity (Cg), 
%         -- Time to travel from x = 0 to actual location (tau)
%         -- Equivalent wavenumber(Keq)
%
%         -- x: this is the space domain. It can be a vector
%         -- d: this is the depth and must be positive. length(d) must be
%         length(x)
%         -- frequency (Hz): Target frequency 

% Author: Enrique M Padilla. 

function [Kinematic] = FUNCTION_Kinematic_FreeWaves(X,h,frequency)

d(1:length(h)) = h; clear h
x(1:length(X)) = X; clear X

if x(1) ~= 0 
    x = horzcat(0,x);
    d = horzcat(d(1),d);
    Kinematic.k = FUNCTION_DispersionEq(d,1/frequency);
    Kinematic.c = (2*pi*frequency)./Kinematic.k; 
    Kinematic.Cg = Kinematic.c./2.*(1+(2.*Kinematic.k.*d)./sinh(2.*Kinematic.k.*d));
    Kinematic.tau = cumsum(abs(horzcat(0,diff(x)))./Kinematic.c);
    Kinematic.Keq = (2*pi*frequency).*Kinematic.tau./abs(x-x(1));
    
    Kinematic.k = Kinematic.k(2:end);
    Kinematic.c = Kinematic.c(2:end);
    Kinematic.Cg = Kinematic.Cg(2:end);
    Kinematic.tau = Kinematic.tau(2:end);
    Kinematic.Keq = Kinematic.Keq(2:end);
else
    Kinematic.k = FUNCTION_DispersionEq(d,1/frequency);
    Kinematic.c = (2*pi*frequency)./Kinematic.k; 
    Kinematic.Cg = Kinematic.c./2.*(1+(2.*Kinematic.k.*d)./sinh(2.*Kinematic.k.*d));
    Kinematic.tau = cumsum(abs(horzcat(0,diff(x)))./Kinematic.c);
    Kinematic.Keq = (2*pi*frequency).*Kinematic.tau./abs(x-x(1));
    Kinematic.Keq(1) = Kinematic.Keq(2)-(Kinematic.Keq(3)-Kinematic.Keq(2));
end
