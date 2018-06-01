function [Ts_mat, name_s] = pretty_print_ts_data(TOL, tol_mode, varargin)
% pretty_print_ts_data(TOL, tol_mode, step_exp1, step_exp2, ...)
  
  name_s = {};
  Ts_mat = [];
  time_scale = 1000; % milliseconds
  name_pad = '    ';
  s_header = '';
  for k=1:length(varargin)
    step_exp_k = varargin{k};
    Ts_vec_k = step_exp_k.settle_time(TOL, tol_mode, 0);
    name_s{k} = step_exp_k.name;
    Ts_mat(:,k) = Ts_vec_k;
    
    s_header = sprintf('%s%s%s', s_header, name_s{k}, name_pad);
  end
  fprintf('%s\n', s_header);
  fprintf('%s\n', repmat('-', 1, length(s_header)-length(name_pad)) );
  
  Ts_mat = Ts_mat*time_scale;

  max_s = floor(max(Ts_mat)); % max of each column.

  for k=1:size(Ts_mat, 1)
    for j = 1:length(name_s)
      s_max = sprintf('%.3f', max_s(j));
      s = sprintf('%.3f', Ts_mat(k,j));
      parts = split(s, '.');
      max_parts = split(s_max, '.');
      int_part = char(parts(1,:));
      int_part_max = char(max_parts(1,:));
      pad_len = length(int_part_max) - length(int_part);
      lpad = repmat(' ', 1, pad_len);
      rpad = repmat(' ', 1, max(length(name_s{j}) - length(s)-length(lpad), 0));
      
      fprintf('%s%s%s%s', lpad, s, rpad, name_pad)
    end
    
    fprintf('\n')
  end
  fprintf('%s\n', repmat('-', 1, length(s_header)-length(name_pad)) );
  
  % Print totals
  for j = 1:length(name_s)
    ts_tot = sum(Ts_mat(:, j));
    
    s = sprintf('%.3f', ts_tot);
    rpad = repmat(' ', 1, length(name_s{j})+length(name_pad) - length(s));
    fprintf('%s%s', s, rpad)
  end
  fprintf('\n')
  
end