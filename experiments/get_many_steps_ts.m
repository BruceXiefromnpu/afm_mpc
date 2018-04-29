
function TS_s = get_many_steps_ts(y, ref_s, step_idx, TOL, verbose, fid)
  % TS_s = get_many_steps_ts(y, ref_s, step_idx, TOL, verbose, fid)
  if ~exist('fid', 'var')
    fid = 1;
  end

  % could do this by reshaping, but in the future, I don't expect all the steps
  % to last for the same amount of time. So we need a cell array.
  Y_cell = {};
  TS_s = zeros(length(ref_s)-1, 1);
  for k = 2:length(ref_s)

    idx_start = step_idx(k);
    if k == length(ref_s)
      idx_end = length(y.Time);
    else
      idx_end = step_idx(k+1)-1;
    end
    idx = idx_start:idx_end;
    Y_cell{k} = timeseries(y.Data(idx), y.Time(idx));

    ref_k = ref_s(k);
    ref_prev = ref_s(k-1);

    delta_ref = ref_k - ref_prev;
    TOL = 0.01*abs(delta_ref);

    if verbose >1
      ts_k = settle_time(Y_cell{k}.Time, Y_cell{k}.Data, ref_k, TOL, k);
    else
      ts_k = settle_time(Y_cell{k}.Time, Y_cell{k}.Data, ref_k, TOL);
    end
    TS_s(k-1) = ts_k;
    if verbose > 0
    fprintf('Ts = %.5f [ms], ref = %.1f, delta_ref = %.1f\n', ts_k*1000,...
            ref_k, delta_ref);
    end
end