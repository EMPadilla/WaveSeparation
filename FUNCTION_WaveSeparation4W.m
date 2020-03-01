%% FUNCTION_WaveSeparation
%  This function separates ingoing and outgoing waves. Thus, this function 
%  may separate incident bound waves, incident free waves and outgoing free 
%  waves along the bathymetry if the propagation chracteristics of the bound
%  waves are properly computed. Note that the configuration of the local
%  array is an imput of the function
%
%             -- signal: It is the time series. It is a matrix where the
%                        rows are the time instants and the columns, the locations.
%             -- X: space domain. This is a vector whose coordinate system
%                   origin is the wave paddles.
%             -- d: This is the water depth and must be positive. length(d) must be
%                   length(X)
%             -- Vbound: This is a vector with the same size as X. Each
%                        element represents the instant celerity of the
%                        bound wave at each X location
%             -- frequency: It is the target frequency where the separation is going
%                           to be performed
%             -- fs: It is the sampling frequency (inverse to Delta t)
%             -- Config: Configuration of the local Array. It is a vector
%                        with 2 values: [P Psep], where P is the 1 sided number
%                        of gauges composing the local array not taking into account 
%                        the reference gauge, i.e., P = (Nx-1)/2, where Nx
%                        is the number of gauges forming the local array.
%                        Conversely. Psep is the jump between gauges with defining 
%                        the resolution of the array, i.e., Psep = 2 means that the 
%                        local array is formed of 1 out of 2 existing gauges. 
%
%             -- IBW: It is an structure containing the following information about the IBW:
%                        - IBW.eta: Time series of IBW
%                        - IBW.Amp: Cross-shore amplitudes along X
%                        - IBW.Ph: Cross-shore phases along X
%                        - IBW.Amp0: Initial Amplitude at X=0 (wave paddle)
%                        - IBW.Ph0: Initial phase at X=0 (wave paddle)
%                        - IBW.alpha: IBW growth rate
%             -- OBW: It is an structure containing the following information about the OBW:
%                        - OBW.eta: Time series of OBW
%                        - OBW.Amp: Cross-shore amplitudes along X
%                        - OBW.Ph: Cross-shore phases along X
%                        - OBW.Amp0: Initial Amplitude at X=0 (wave paddle)
%                        - OBW.Ph0: Initial phase at X=0 (wave paddle)
%                        - OBW.alpha: OBW growth rate
%             -- OFW: It is an structure containing the following information about the OFW:
%                        - OFW.eta: Time series of OFW
%                        - OFW.Amp: Cross-shore amplitudes along X
%                        - OFW.Ph: Cross-shore phases along X
%                        - OFW.Amp0: Initial Amplitude at X=0 (wave paddle)
%                        - OFW.Ph0: Initial phase at X=0 (wave paddle)
%             -- IFW: It is an structure containing the following information about the IFW:
%                        - IFW.eta: Time series of IFW
%                        - IFW.Amp: Cross-shore amplitudes along X
%                        - IFW.Ph: Cross-shore phases along X
%                        - IFW.Amp0: Initial Amplitude at X=0 (wave paddle)
%                        - IFW.Ph0: Initial phase at X=0 (wave paddle)
%
% Author: Enrique M. Padilla

function [IBW,OBW,OFW,IFW] = FUNCTION_WaveSeparation4W(signal,X,d,VBound,frequency,fs,Config)
%% Definition of the array
Nx = Config(1);
Psep = Config(2); 
Array = 1:Psep:Nx*Psep-(Psep-1);
ref = ceil(Nx/2);

%% Defines the shoreline limit, since the separation is not valid beyond the shoreline

if d(1) * d(end) <= 0 % If this number is negative, there a sholine somewhere in th domain
    pshore = FUNCTION_Closest(0,d);
else 
    pshore = length(X);
end
shoreline = X(pshore);

%% Propagation parameters for the Bound Wave

% needs propagation from X=0
if X(1)==0
    Kinematic.BoundAprox.tau = cumsum(abs(horzcat(0,diff(X)))./VBound);
    Kinematic.BoundAprox.Keq = (2*pi*frequency).*Kinematic.BoundAprox.tau./abs(X-X(1));
    Kinematic.BoundAprox.Keq(1) = Kinematic.BoundAprox.Keq(2);
else
    Kinematic.BoundAprox.tau = cumsum(abs(horzcat(0,diff(horzcat(0,X))))./horzcat(VBound(1),VBound));
    Kinematic.BoundAprox.Keq = (2*pi*frequency).*Kinematic.BoundAprox.tau./abs(horzcat(0,X));

    Kinematic.BoundAprox.tau = Kinematic.BoundAprox.tau(2:end);
    Kinematic.BoundAprox.Keq = Kinematic.BoundAprox.Keq(2:end);    
end

%% Kinematics for free waves (F = fg). Linear Theory
[Kinematic.Free] = FUNCTION_Kinematic_FreeWaves(X,d,frequency);

%% frequency window around the target frequency 
[Pxx,~,F] = FUNCTION_Fourier_IntegralSpectrum(signal(:,1),fs);
frequency = F(abs((F-frequency))==min(abs(F-frequency)));
Boundary = FUNCTION_PeakThreshold(Pxx(:,1),frequency,F);
window = F(Boundary);

%% Estimation of the complex coefficients of the signals
Lt = length(signal(:,1));
LX = length(X);

Gn = fft(signal,Lt,1);
Gn(1,:) = 0;

Z_IBW = zeros(Lt,LX);
Z_OFW = zeros(Lt,LX);
Z_IFW = zeros(Lt,LX);
Z_OBW = zeros(Lt,LX);

%% Time domain
vart = 1/fs;
t = 0:vart:vart*Lt-vart;

%% Separation 
alphaIBW(1) = 1;
alphaOBW(1) = 1;
for It = 1:5
    for loc = 1:pshore
        vec = loc-Array(ref)+Array;
        vec = vec(vec > 0 & vec < pshore);
        while length(vec) < 4
            vec = horzcat(vec(1)-Psep,vec,vec(end)+Psep);
            vec = vec(vec > 0 & vec < pshore);
        end       
        % Propagation coefficients               
        Q_IBW = (d(loc)./d(vec)).^(alphaIBW(It)) .* ...
                exp(-1i*(Kinematic.BoundAprox.Keq(vec).*(X(vec))-Kinematic.BoundAprox.Keq(loc)*(X(loc))));       
        Q_OBW = (d(loc)./d(vec)).^(alphaOBW(It)) .* ...
                exp(1i*(Kinematic.BoundAprox.Keq(vec).*(X(vec))-Kinematic.BoundAprox.Keq(loc)*(X(loc))));                 
        Q_IFW = sqrt(Kinematic.Free.Cg(loc)./Kinematic.Free.Cg(vec)) .* ...
                exp(-1i*(Kinematic.Free.Keq(vec).*(X(vec))-Kinematic.Free.Keq(loc)*(X(loc))));
 
        Q_OFW  = sqrt(Kinematic.Free.Cg(loc)./Kinematic.Free.Cg(vec)) .* ...
                exp(1i*(Kinematic.Free.Keq(vec).*(X(vec))-Kinematic.Free.Keq(loc)*(X(loc))));

        %Solving the system
        A(:,1) = Q_IBW;
        A(:,2) = Q_OFW;
        A(:,3) = Q_IFW;
        A(:,4) = Q_OBW;
        for nf = find(F==window(1)):find(F==window(2))
            b(:,1) = Gn(nf,vec);

            Z = A \ b;
            Z_IBW(nf,loc) = Z(1);
            Z_OFW(nf,loc) = Z(2);
            Z_IFW(nf,loc) = Z(3);
            Z_OBW(nf,loc) = Z(4);
            
            Z_IBW(Lt+2-nf,loc) = conj(Z_IBW(nf,loc));
            Z_OFW(Lt+2-nf,loc) = conj(Z_OFW(nf,loc));
            Z_IFW(Lt+2-nf,loc) = conj(Z_IFW(nf,loc));
            Z_OBW(Lt+2-nf,loc) = conj(Z_OBW(nf,loc));
        end
        clear A b
    end

    % Build the time series from the complex amplitudes
    IBW.eta = real(ifft(Z_IBW,Lt,1));
    OFW.eta = real(ifft(Z_OFW,Lt,1)); 
    IFW.eta = real(ifft(Z_IFW,Lt,1));
    OBW.eta = real(ifft(Z_OBW,Lt,1));


    % Iterate the alpha value for IBW
    [IBW.Amp,~] = FUNCTION_AmplitudePhase(IBW.eta,t,frequency); close
    BoundSh_functionIBW = @(coef,xdata) IBW.Amp(1)*(d(1)./xdata).^(coef(1));
    CC = lsqcurvefit(BoundSh_functionIBW, alphaIBW(It), d, IBW.Amp);
    alphaIBW(It+1) = CC(1);
    
    % Iterate the alpha value for OBW
    [OBW.Amp,~] = FUNCTION_AmplitudePhase(OBW.eta,t,frequency); close
    BoundSh_functionOBW = @(coef,xdata) OBW.Amp(1)*(d(1)./xdata).^(coef(1));
    DD = lsqcurvefit(BoundSh_functionOBW, alphaOBW(It), d, OBW.Amp);
    alphaOBW(It+1) = DD(1);
end


%% Outcomes

% alpha and beta values for ILW
IBW.alpha = alphaIBW;
OBW.alpha = alphaOBW;
%plot(X, IBW.Amp,'k.');hold on;plot(X,BoundSh_functionIBW(alphaIBW(It+1),d),'r-')
%plot(X, OBW.Amp,'b.');hold on;plot(X,BoundSh_functionOBW(alphaOBW(It+1),d),'g-')

% Amplitudes and Phases
[IBW.Amp,IBW.Ph] = FUNCTION_AmplitudePhase(IBW.eta,t,frequency); close
[OFW.Amp,OFW.Ph] = FUNCTION_AmplitudePhase(OFW.eta,t,frequency); close
[IFW.Amp,IFW.Ph] = FUNCTION_AmplitudePhase(IFW.eta,t,frequency); close
[OBW.Amp,OBW.Ph] = FUNCTION_AmplitudePhase(OBW.eta,t,frequency); close

% Initial Amplitudes and Phases
IBW.Amp0 = IBW.Amp(1);
IBW.Ph0 = mod(IBW.Ph(ref)+Kinematic.BoundAprox.Keq(ref)*X(ref),2*pi);

OBW.Amp0 = OBW.Amp(1);
OBW.Ph0 = mod(OBW.Ph(ref)-Kinematic.BoundAprox.Keq(ref)*X(ref),2*pi);

IFW.Amp0 = IFW.Amp(ref).* sqrt(Kinematic.Free.Cg(ref)./Kinematic.Free.Cg(1));
IFW.Ph0 = mod(IFW.Ph(ref)+Kinematic.Free.Keq(ref)*X(ref),2*pi);

OFW.Amp0 = OFW.Amp(ref).* sqrt(Kinematic.Free.Cg(ref)./Kinematic.Free.Cg(1));
OFW.Ph0 = mod(OFW.Ph(ref)-Kinematic.Free.Keq(ref)*X(ref),2*pi);






function [Boundary] = FUNCTION_PeakThreshold(signal,Peak,F)
Blow = 1;
Bhigh = length(signal);
PositionPeak = find(abs((F-Peak))==min(abs(F-Peak)));
lim = 0.1; % 0.1 Hz is the limit

cont = 1;
while Blow == 1 && PositionPeak-cont > 0
    if cont==PositionPeak
       Blow = 1;
    else
        if (signal(PositionPeak-cont+1)-signal(PositionPeak-cont))/signal(PositionPeak-cont+1)<=0.00001 || abs(F(PositionPeak-cont)-Peak)>lim
            Blow = PositionPeak-(cont-1);
        else
            cont = cont+1;
        end
    end
end
     
cont = 1;
while Bhigh==length(signal) && PositionPeak+cont < length(signal)
    if (signal(PositionPeak+cont-1)-signal(PositionPeak+cont))/signal(PositionPeak+cont-1)<=0.00001 || abs(F(PositionPeak+cont)-Peak)>lim
        Bhigh = PositionPeak+(cont-1);
    else
        cont = cont+1;
    end
end

Boundary(1) = Blow;
Boundary(2) = Bhigh;



function [pos] = FUNCTION_Closest(vector1,vector2) % for each value in vector 1, look for the closest in vector2
    for n = 1:length(vector1)
        pos(n) = min(find(abs((vector2-vector1(n)))==min(abs(vector2-vector1(n)))));
    end