%% FUNCTION_DispersionEq
%  This function provides the wavenumber solving the dispersion equation
%  for a given period and water depth
%
%         -- h: this is the depth and must be positive. h can be a vector
%               for different locations if Period is not a vector
%         -- Period (s): is the inverse of the target frequency in Hz. Period can 
%                        be a vector in case multiple frequencies want to be computed 
%                        if h is not a vector
%         -- k: is the wavenumber. 

% Author: Enrique M Padilla. 

function [k] = FUNCTION_DispersionEq(h,Period)

accuracy = 0.01;
lambda(:,1) = 0:accuracy:1000;
if length(h) > 1
   Period = ones(1,length(h))*Period;
else
   h = ones(1,length(Period))*h; 
end
eq = abs((lambda/(2*pi*9.81))*(2*pi./Period).^2-tanh(2*pi*1./lambda*h));
[Eq,Lambda] = meshgrid(min(eq,[],1),lambda);
k(1,:) = 2*pi./Lambda(eq == Eq);
