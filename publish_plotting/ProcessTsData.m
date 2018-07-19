classdef ProcessTsData

  methods (Static)
    
    function S = TS_dat2tex(TS_dat_cell, step_ref, varargin)
      % S = TS_dat2tex(TS_dat_cell, step_ref, varargin)
      % S = TS_dat2tex(TS_dat_cell, step_ref, 'do_color', (true|false))
      % S = TS_dat2tex(TS_dat_cell, step_ref, 'do_color', (true|false), ts_vec)
      %
      % Purpose:
      % -------
      %  To format settling time data from a bunch of step experiments into a latex
      %  table. This table will be formatted as each row corresponding to a different
      %  step input, and each column corresponding to a different experiment (e.g.,
      %  col1 is linear feedback, col2 is MPC).
      %
      % Inputs
      % ------
      %   TS_dat_cell : A cell array containing structures with the expected fields
      %                 of TS_dat_cell{i}.ts_s and TS_dat_dat_cell{i}.name
      %                 where ts_s is a vector of settling times for the ith
      %                 experiment, name is a string containing the name of the
      %                 experiment. The vector ts_s should be in the same order as
      %                 the steps in step_ref
      %
      %   step_ref : an instance of the class StepRef which should contain the step
      %              reference data for the experiments.
      %
      % Optional inputs (name-value pairs)
      %  --------------
      %    'do_color', (true|false) : Wheather or not to color the background of the
      %    latex table as a color map corresponding to the length of the settle time.
      %    'ts_vec', ts_vec : a vector of total settling times, to be used in
      %    computing the color map. This is useful if you are building another table
      %    with different data and need both colormaps to have the same range. There
      %    is probably a better way but this works for now.
      %
      %  Outputs
      %  -------
      %  S : A (large!) string containing the the latex table. Use fprintf to write
      %  this to a .tex file.
      
      p = inputParser();
      p.addParameter('do_color', true);
      p.addParameter('ts_vec', []);
      
      p.parse(varargin{:});
      do_color = p.Results.do_color;
      ts_vec   = p.Results.ts_vec;
      
      % ------------------------------------------------------
      % ------------ Build the LaTex table -------------------
      
      if do_color
        % Build a colormap, interpolated onto the TOTAL number of settling times we have.
        mp = colormap(gca, 'jet');
        x = 1:size(mp,1);
        if isempty(ts_vec)
          ts_vec = concat_TS_cell(TS_dat_cell);
        end
        ts_vec = unique(sort(ts_vec, 'descend'));
        xq = linspace(1, (size(mp,1)), length(ts_vec));
        
        mp_fine = zeros(length(xq), 3);
        mp_fine(:,1) = interp1(x, mp(:,1), xq);
        mp_fine(:,2) = interp1(x, mp(:,2), xq);
        mp_fine(:,3) = interp1(x, mp(:,3), xq);
        
      end
      
      % -- First, we programmatically construct \tabular{ccc},
      %    since the 'ccc' depends on how many columns we need.
      %    the number of columns in the table
      c_fmt = repmat('c', 1, 3+length(TS_dat_cell));
      S = sprintf('\\begin{tabular}{%s}\n', c_fmt);
      
      % -- Form the table header:
      str_ref_cols = sprintf('&ref & delta');
      str_dat_cols = '';
      for k=1:length(TS_dat_cell)
        str_dat_cols = sprintf(' %s & %s', str_dat_cols, TS_dat_cell{k}.name);
      end
      S = sprintf('%s%s%s\\\\\n\\toprule\n', S, str_ref_cols, str_dat_cols);
      %
      % -- Build up the body of the table. Outer loop is for each row.
      % Inner loop is for each experiment (columns).
      
      
      for k = 2:length(step_ref.step_amps)
        
        delta_ref = step_ref.step_amps(k) - step_ref.step_amps(k-1);
        str_ref_cols = sprintf('&%.2f & %.2f', step_ref.step_amps(k), delta_ref);
        str_dat_cols = '';
        % across columns, each experiment
        for j = 1:length(TS_dat_cell)
          % settle times are stored as columns, but we have to build each row.
          ts_kj = TS_dat_cell{j}.ts_s(k-1); % kth row, jth col.
          
          % Indexes from slowest to fastest.
          if do_color
            idx = find(ts_vec == ts_kj);
            str_dat_cols = sprintf('%s &\\cellcolor[rgb]{%.4f, %.4f, %.4f} %.2f', str_dat_cols, ...
              mp_fine(idx, 1), mp_fine(idx,2), mp_fine(idx, 3), 1000*ts_kj);
           % keyboard
          else
            str_dat_cols = sprintf('%s & %.2f', str_dat_cols,  1000*ts_kj);
          end
          % keyboard
        end
        s_row = sprintf('%s%s\\\\ \n', str_ref_cols, str_dat_cols);
        S = sprintf('%s%s', S, s_row);
      end
      % -- Build the footer. This is where the totals go.
      str_ref_cols = sprintf('total & -- & --');
      str_dat_cols = '';
      for j = 1:length(TS_dat_cell)
        str_dat_cols = sprintf('%s &%.2f', str_dat_cols, ...
          1000*sum(TS_dat_cell{j}.ts_s) );
      end
      
      s_row = sprintf('%s%s\\\\ \n', str_ref_cols, str_dat_cols);
      
      S = sprintf('%s\\midrule\n %s', S, s_row);
      
      
      S = sprintf('%s\\end{tabular}\n', S);
      
    end
    function ts_vec = ts_vec_from_dir(root, TOL, tol_mode)
      files = strsplit(ls(root));
      
      ts_vec = [];
      verbose = 0;
      for file = files
        file_path = fullfile(root, file{1});
        if ~isfile(file_path)
          continue; % skip directories.
        end
        dat = load(file_path);
        exp_name_str = fields(dat);
        exp_name_str = exp_name_str{1};
        dat = dat.(exp_name_str);
        ts_vec = [ts_vec; dat.settle_time(TOL, tol_mode, verbose)];
      end
    end
    function ts_vec = concat_TS_cell(ts_cell)
      ts_vec = []
      for k=1:length(ts_cell)
        ts_vec = [ts_vec; ts_cell{k}.ts_s];
      end
    end
    
  end
end