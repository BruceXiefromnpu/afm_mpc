function [ y, y_vec_k ] = hyst_play_op(u, r, w_y0)
  %HYST_PLAY_OP Summary of this function goes here
  %   Detailed explanation goes here
  
  
  n = length(r);
  w = w_y0(1:n);
%   y0 = w_y0(n+1:end);
  y0 = w*0;
  
  y_vec_k = zeros(length(u), length(r));
  y = 0*u;

  y_vec_k(1, :) = y0(:)';

  for k=2:length(u)
    uk = u(k);
    for j = 1:length(r)
      y_vec_k(k, j) = max(uk - r(j), min(uk+r(j), y_vec_k(k-1, j)));
    end
    y(k) = w'*y_vec_k(k,:)';
  end
  
  
  
end

