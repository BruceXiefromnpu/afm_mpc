function y_out = rate_limit(y, dy_max)
% y_out = rate_limit(y, dy_max)
% Enforces a rate of change limit of dy_max on the input vector y.
  y_out = y*0;
  y_out(1) = y(1);
  
  for k = 2:length(y)
    
    dy = y(k) - y_out(k-1);
    
    if dy > 0
      dy = min(dy_max, dy);
    else
      dy = max(-dy_max, dy);
    end
    y_out(k) = y_out(k-1) + dy;
    
  end
  
  
end