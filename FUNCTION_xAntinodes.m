%% FUNCTION_xAntinodes
%  This function provides the expected location of nodes and antinodes resultant
%  of 2 pair of wave trains travelling with different celerities at the
%  same frequency
%
%         -- x: Spatial domain. It is a vector
%         -- f: Target frequency in Hz 
%         -- phases: Vector containing the initial pahses (in radians) for the 2 wave trains, e.g., phases = [0.3 3.4] 
%         -- celerities: celerity vector for each wave train along the
%                        spatial domain. Therefore celerities(1,:) is the
%                        celerity vector for the wave train 1. Idem for
%                        wave train 2
%
%         -- X_anti: Cross-shore locations where the antinodes are expected
%         -- X_nodes: Cross-shore locations where the nodes are expected
%         -- L_Approx: Flat bed approximation of the distance between
%                      consecutive nodes

% Author: Enrique M Padilla.


function [X_anti,X_nodes,L_Approx] = FUNCTION_xAntinodes(x,f,phases,celerities)
c1 = celerities(1,:);
c2 = celerities(2,:);

% x must begin in x = 0;
if x(1) ~= 0
    x = horzcat(0,x);
    c1 = horzcat(c1(1),c1);
    c2 = horzcat(c2(1),c2);
end

%% Interpolate vectors
dX = 0.01;
Xinterp = 0:dX:x(end);
C1interp = interp1(x,c1,Xinterp);
C2interp = interp1(x,c2,Xinterp);

%% Searching for nodes and antinodes

for nx = 2:length(Xinterp)
    int(nx) = (2*pi*f)*(trapz(Xinterp(1:nx),(1./C1interp(1:nx)-1./C2interp(1:nx))));
end

m0 = -100; % this tries only -100<m<100
X_anti = [];
X_nodes = [];
for m = m0:-m0
    function_antinodes = int+(phases(1)-phases(2))+(2*pi*m);
    sg = min(function_antinodes(1) * function_antinodes);
    if sg < 0     
        X_anti = horzcat(X_anti,Xinterp(find(abs((function_antinodes-0))==min(abs(function_antinodes-0)), 1 )));
    end
    function_nodes = int+(phases(1)-phases(2)) + pi*(2*m+1);
    sg = min(function_nodes(1) * function_nodes);
    if sg < 0     
        X_nodes = horzcat(X_nodes,Xinterp(find(abs((function_nodes-0))==min(abs(function_nodes-0)), 1 )));
    end
end

%% Approximation for flat bed
L_Approx = (1/f)*abs((c1(1)*c2(1))/(c2(1)-c1(1)));

    
