function [pow_spec, freqs] = power_spectrum(x, Ts )
  % [pow_spec, freqs] = power_spectrum(x, Ts )
  % For a vector of data x at sample rate Ts,
  % returns the power spectrum and associated frequencies.
  % Roughly this is the modulus of the first half of the FFT of the signal.

  Fs = 1/Ts;

  Y = fft(x);
  L = length(x);
  P2 = abs(Y/L);
  pow_spec = P2(1:L/2+1);
  pow_spec(2:end-1) = 2*pow_spec(2:end-1);
  freqs = Fs*[0:L/2]/L;

end

