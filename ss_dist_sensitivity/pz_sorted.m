function [p_sorted, z_sorted] = pz_sorted(sys)
  if sys.Ts == 0
    error('This function is only implemented for discrete time systems')
  end
  p = pole(sys);
  z = tzero(sys);
  
  % sort pole by frequency:
  sp_locs = log(p)/sys.Ts;
  p_wns = abs(sp_locs);
  [~, idx_p] = sort(p_wns);
end

