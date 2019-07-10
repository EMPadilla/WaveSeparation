%% FUNCTION_AmplitudePhase
%  This function provides the amplitude and phase of a time
%  series for a range of locations at a certain frequency.
%
%           -- Signal: is the time series. It may be matrix where every row
%                      is the time series at different locations (each column)
%           -- Time: time domain (vector)
%           -- frequency: this is the target frequency (in Hz) where the amplitude and
%                         phase are computed.
%           -- Amplitude: Vector of Amplitudes. Each element is the amplitude at every location.
%           -- Phase: Vector of Phases. Each element is the Phase at every location.
%
% Author: Enrique M Padilla. 

function [Amplitude,Phase] = FUNCTION_AmplitudePhase(Signal,Time,frequency) 
Lx = length(Signal(1,:));
vart = Time(2)-Time(1);
fs = round(1/vart);

[Pxx,~,F] = FUNCTION_Fourier_IntegralSpectrum(Signal,fs);  % Wave Spectra (Pxx) and frequency domain (F)
[~,loc] = findpeaks(Pxx(:,1));close % locate the positions in the freqeuncy domain with energetic peaks 
if isempty(loc) % if no peaks (when signal is a zero-vector)
    Amplitude = zeros(1,Lx);
    Phase = zeros(1,Lx);
    sprintf('Energy at frequency %s Hz is null. Check if expected',num2str(frequency))  
else
    Floc = F(loc);
    frequency = Floc(find(abs((Floc-frequency))==min(abs(Floc-frequency)),1)); % closer frequency value within the frequency domain
    
    Boundary = FUNCTION_PeakThreshold(Pxx(:,1),frequency,F); % boundary around the target frequency
    if Boundary(1) == Boundary(2)
        Boundary(1) = Boundary(1)-1;
        Boundary(2) = Boundary(2)+1;
    end
    Amplitude = sqrt(2*trapz(F(Boundary(1):Boundary(2)),Pxx(Boundary(1):Boundary(2),:)));
    signal_filt = FUNCTION_Filtering_Components([F(Boundary(1)) F(Boundary(2))],Signal,fs);close
    Phase = FUNCTION_Phase(signal_filt,Time,frequency);

    Energy = trapz(F(Boundary(1):Boundary(2)),Pxx(Boundary(1):Boundary(2),:));
    EnergyTotal = trapz(F(2:end),Pxx(2:end,:));
    sprintf('Energy at frequency %s Hz comprises a %s %% of the total energy content',num2str(frequency),num2str(Energy*100/EnergyTotal))   
end

%% Checking graph 
graph = 1;
if graph == 1 && ~isempty(loc)
   subplot(2,1,1)
      plot(Time,signal_filt(:,1),'.');hold on
      plot(Time,Amplitude(1)*cos(2*pi*frequency*Time+Phase(1)),'r')
      plot(Time(round(length(Time)/2)),signal_filt(round(length(Time)/2),1),'ro')
        
   subplot(2,1,2)
       plot(F(2:end),Pxx(2:end,1));hold on
       plot([F(Boundary(1)) F(Boundary(2))],[Pxx(Boundary(1),1) Pxx(Boundary(2),1)],'ro')
       xlim([F(Boundary(1))*0.95 F(Boundary(2))*1.05])
       xlabel('Frequencies (Hz)')
end




% Functions ****************************************************
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

if Blow < PositionPeak-10
   Blow = PositionPeak-10;
end

if Bhigh > PositionPeak+10
   Bhigh = PositionPeak+10;
end

Boundary(1) = Blow;
Boundary(2) = Bhigh;

function [Phase] = FUNCTION_Phase(signal,Time,frequency)

Phase = ones(size(signal(1,:)))*NaN;
vart = Time(2)-Time(1);

[m,~] = meshgrid(mean(signal,1),1:length(signal(:,1)));
[c,~] = meshgrid((1./max(signal-m,[],1)),1:length(signal(:,1)));
Signal = (signal-m).*c;

N = round(length(signal)/2);
for loc = 1:length(signal(1,:))
    if Signal(N,loc)-Signal(N+1,loc) > 0
       PhaseN = acos(Signal(N,loc));
       Phase(loc) = PhaseN - 2*pi*frequency*(mod(N*vart,1/frequency));
    end
    if Signal(N,loc)-Signal(N+1,loc) < 0
       PhaseN = 2*pi-acos(Signal(N,loc));
       Phase(loc) = PhaseN - 2*pi*frequency*(mod(N*vart,1/frequency));
    end
    if Signal(N,loc)-Signal(N+1,loc) == 0
       Phase(loc) = 0;
    end
end

function [Spectrum,Phase,F] = FUNCTION_Fourier_IntegralSpectrum(signal,fs)
N = length(signal(:,1));
Z = fft(signal,N,1);

if mod(N,2) == 0
    G = Z(1:N/2,:);
else
    G = Z(1:(N+1)/2,:);
end

Spectrum = (1/(fs*N)).*abs(G).^2;
Spectrum(2:end-1,:) = 2*Spectrum(2:end-1,:);
F = linspace(0,fs/2,length(G(:,1)));

Phase = unwrap(angle(G),[],1);
