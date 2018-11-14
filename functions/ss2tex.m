function str = ss2tex(sys)



  g = zpk(sys);
  [z, p, k] = zpkdata(sys, 'v');
  
str = sprintf('\\begin{tabular}{cc}\n');
str = sprintf('%spole & zero\\\\\n', str);

p_real = p(imag(p)==0);
z_real = z(imag(z)==0);

p_imag = zp_noconj(p);
z_imag = zp_noconj(z);

P_imag = count_mults(p_imag);
Z_imag = count_mults(z_imag);

P_real = count_mults(p_real);
Z_real = count_mults(z_real);

for k = 1:length(P_imag)
  if k <=size(Z_imag,1)
    str = sprintf('%s%.4g$\\pm$%.4f & %.4f $\\pm$ %.4f\\\\ \n', str,...
      real(P_imag(k,1)), imag(P_imag(k,1)), real(Z_imag(k,1)), imag(Z_imag(k,1)));
  else
    str = sprintf('%s%.4g$\\pm$%.4f & --\\\\ \n', str,real(P_imag(k,1)), imag(P_imag(k,1)));
  end
  
end
  
for k = 1:length(P_real)
  p_str = real_mult_str(P_real, k);
  
  if k <=size(Z_real,1)
    z_str = real_mult_str(Z_real, k);
  else
    z_str = '--';
  end
    str = sprintf('%s%s & %s\\\\ \n', str, p_str, z_str);
end


str = sprintf('%s\\end{tabular}', str);
end


function str = real_mult_str(PZ_real, row_idx)
% PZ_real should be a matrix containing the pole OR zero locations in the first
% column, and their multiplicities in the second column.
  if PZ_real(row_idx,2) == 1
    str = sprintf('%.4g', PZ_real(row_idx, 1));
  elseif PZ_real(row_idx,2) > 1
    str = sprintf('(%d) %.4g', PZ_real(row_idx, 2), PZ_real(row_idx, 1));
  else
    str = ''
  end
    
end



function zp_imag = zp_noconj(zp)
zp_imag = [];

  while ~isempty(zp)
    if imag(zp(1)) == 0
      zp(1) = [];
      continue
    end

    zp_imag = [zp_imag; zp(1)];
    zp(1) = [];

    [p_conj, idx] = find(zp == conj(zp_imag(end)), 1, 'first');
    zp(idx) = [];

  end
  
  
end

function ZP_mults = count_mults(zp)
  ZP_mults = [];

  while ~isempty(zp)
    
    p_i = zp(1);
    p_all = zp(zp == p_i);
    zp( zp == p_i) = [];
    ZP_mults = [ZP_mults; p_i, length(p_all)];
    
  end
  
%   P_mults

end