function [P1, freqs] = power_spectrum(x, Ts )
%FFT_SPEC Summary of this function goes here
%   Detailed explanation goes here

Fs = 1/Ts;

Y = fft(x);
L = length(x);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
freqs = Fs*[0:L/2]/L;

end

