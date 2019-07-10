%% FUNCTION_Filtering_Components
%  This function provides a filtered signal for a desired frequency
%
%             -- window: It is a vector of 2 values: the lower and upper 
%                        frequency limit where the signal is going to be filter. 
%             -- signal: It is the time series. It may be a vector where the
%                        rows are the time series and the columns, the locations. 
%             -- fs: is the sampling frequency. This is inverse to Delta t
%             -- signal_filt: is the filtered time series at the frequency
%                             defined by the window. The format is the same
%                             as signal

% Author: Enrique M. Padilla

function [signal_filt] = FUNCTION_Filtering_Components(window,signal,fs)

Lt = length(signal(:,1));
if mod(Lt,2) == 0
    F = linspace(0,fs/2,Lt/2);
else
    F = linspace(0,fs/2,(Lt+1)/2);
end

Gn = fft(signal,Lt,1);
LGn = length(Gn(:,1));

F_low = window(1);
F_high = window(2);

Gn([find(F<F_low) find(F>F_high) LGn+2-find(F<F_low)  LGn+2-find(F>F_high)],:) = 0;
signal_filt = real(ifft(Gn,Lt,1));

%% Graphs
Graphs = 1;
if Graphs == 1
    plot(signal(:,1));hold on
    plot(signal_filt(:,1),'.r')
    legend('original signal','filtered signal')
end

