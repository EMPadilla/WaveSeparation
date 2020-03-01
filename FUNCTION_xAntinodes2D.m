%% FUNCTION_xAntinodes2D
%  This function provides the expected location of nodes and antinodes resultant
%  of 2 pair of wave trains travelling with different celerities at the
%  same frequency
%
%         -- x: Spatial domain. It is a vector 
%         -- phases: Vector containing the initial pahses (in radians) for the 2 wave trains, e.g., phases = [0.3 3.4] 
%         -- WaveNumbers: Wave Number vector for each wave train along the
%                        spatial domain. Therefore WaveNumbers(1,:) is the
%                        wave number vector for the wave train 1. Idem for
%                        wave train 2
%
%         -- Anti: Cross-shore locations where the antinodes are expected
%         -- Nodes: Cross-shore locations where the nodes are expected
%         -- L_Approx: Flat bed approximation of the distance between
%                      consecutive nodes

% Author: Enrique M Padilla.

function [Anti,Nodes,L_Approx] = FUNCTION_xAntinodes2D(x,phases,WaveNumbers)

%% interpolate vectors
dX = 0.01;
Xinterp = 0:dX:x(end);
K1 = interp1(x,WaveNumbers(1,:),Xinterp);
K2 = interp1(x,WaveNumbers(2,:),Xinterp);

Anti = [];
Nodes = [];
m0 = -100; % this tries only -100<m<100
for m = m0:-m0 
    function_antinodes = (2*pi*m-(phases(1)-phases(2)))./-(K1(1,:)-K2(1,:))-Xinterp;
    sg = function_antinodes(1) * function_antinodes(end);
    if sg < 0     
        Anti = horzcat(Anti,Xinterp(find(abs((function_antinodes))==min(abs(function_antinodes)), 1 )));
    end
    
    function_nodes = (pi*(2*m+1)-(phases(1)-phases(2)))./-(K1(1,:)-K2(1,:))-Xinterp;
    sg = function_nodes(1) * function_nodes(end);
    if sg < 0     
        Nodes = horzcat(Nodes,Xinterp(find(abs((function_nodes))==min(abs(function_nodes)), 1 )));
    end
end

Anti = sort(Anti);
Nodes = sort(Nodes);

L_Approx = abs(-2*pi./(WaveNumbers(1,1)-WaveNumbers(2,1)));

    
