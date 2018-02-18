% [yy_imp ] = frf2impulse(Gz_even, N)
%
% Calculates impulse response from frf data, Gz_even.
% The function assumes that Gz_even is the FRF at evenly spaced frequency
% from to nyquist.
%
%
% Example: 
%  freqs_even = [0:(1/(N*Ts)):floor((1/Ts)*0.5)]'; %Frequencies up to nyquist.
%  ws_even    = freqs_even*2*pi;
%  K_freq_max = find(ws_even <= w_s(end), 1, 'last')  
%  Gz_even = interp1(w_s, Gz_frf, ws_even(1:K_freq_max), 'spline');
% 
% Then (after delay removed if needed, can do:
% Gz_even(K_freq_max:length(ws_even)) = 0;
%
% 

function [yy_imp ] = frf2impulse(Gz_even)

% yyFF = ([Gz_even;  zeros(length(Gz_even), 1)]);
Gz_even_conj = conj(flipud(Gz_even));

% repeated zero frequency, go to end-1;
yyFF = ([Gz_even;  Gz_even_conj(1:end-1)])/2;

yy_expIFT  = ifft(yyFF);
yy_imp     = 2*real(yy_expIFT);

end

