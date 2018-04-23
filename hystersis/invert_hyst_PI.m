


function [rp, wp] = invert_hyst_PI(r, w)
  r = r(:);
  w = w(:);
  rp = r*0;
  wp = w*0;
  
  wp(1) = 1/w(1);
  
  for i=2:length(w)
    
    s1 = sum(w(1:i));
    s2 = sum(w(1:i-1));
    
    wp(i) = -w(i)/(s1*s2);
    
  end
  
  for i=1:length(r)
    %keyboard
    rp(i) = sum( w(1:i).*(r(i)-r(1:i)));
%     sum(wp(1:i))*r(i) - sum( wp(1:i).*r(1:i));
    
  end
  
  
  
end