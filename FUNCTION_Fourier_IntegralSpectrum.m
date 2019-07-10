%% FUNCTION_Fourier_IntegralSpectrum
%  This function provides the power spectrum of a signal or a group of
%  signals (signal may be  matrix).
%
%          -- signal: It is the time series. It may be a vector where the
%                     rows are the time series and the columns, the locations. 
%          -- fs: sampling frequency. Inverse of Delta t
%          -- Spectrum: One side Energy density spectrum. Units (m^2/Hz)
%          -- Angle: this is the initial phase for every frequency of the
%                    vector F
%          -- F: one side frequency vector (Hz)

% Author: Enrique M Padilla. 

function [Spectrum,Phase,F] = FUNCTION_Fourier_IntegralSpectrum(signal,fs,varargin)
if nargin > 2 % if df needs imposing
   df = varargin{1};
   N = round(fs/df);
else
   N = length(signal(:,1));
end

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

%% Note this is equal to 
% [P Fr] = periodogram(signal,rectwin(length(signal)),length(signal),Fs);
% and equal to
% [P Fr] = pwelch(signal,ones(length(signal),1),[],length(signal),fs);
 
%plot(F,Spectrum);
%hold on; plot(F,Spectrum,'o');plot(Fr,P,'r');plot(Fr,P,'r.')

