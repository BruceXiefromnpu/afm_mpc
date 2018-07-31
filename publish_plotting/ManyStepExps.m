classdef ManyStepExps < handle
  properties
    TOL;
    tol_mode;
    step_ref;
    step_exps;
    TS_mat;
  end
  methods
    function self = ManyStepExps(TOL, tol_mode, step_ref, varargin)
    % ProcessTsData(TOL, tol_mode, step_ref, varargin)
    % varargin should be:
    %  step_exp1, step_exp2,...
      self.TOL = TOL;
      self.tol_mode = tol_mode;
      self.step_ref = step_ref;
      
      for k=1:length(varargin)
        self.step_exps{k} = varargin{k};
        Ts_vec_k = self.step_exps{k}.settle_time(TOL, tol_mode, 0);
        self.TS_mat(:,k) = Ts_vec_k;
      end
    end

    function self = add_step_exp(self, step_exp)
    % self = add_step_exp(self, step_exp)
    % Add another step experiment to the property (a cell array) 
    % self.step_exps
      n = length(self.step_exps);
      self.step_exps(n+1) = step_exp;
    end
    
    function [hands] = ploty_all(self, fig_ax)
    % [hands] = ploty_all(self, ax)
    % plot all of the step experiment y (output) trajectories in
    % self.step_exps to the figure or axis pointed to by fig_ax.
      idx = 1:length(self.step_exps);
      hands = self.ploty_selected(idx, fig_ax);
    end
    function [hands] = ploty_selected(self, idx, fig_ax)
      hands = gobjects(1, length(idx));
      for k=1:length(idx)
        hands(k) = self.step_exps{idx(k)}.ploty(fig_ax);
        hold on
      end
    end
    
    function [hands] = plotdu_selected(self, idx, fig_ax)
      hands = gobjects(1, length(idx));
      for k=1:length(idx)
        hands(k) = self.step_exps{idx(k)}.plotdu(fig_ax);
        hold on
      end
    end
    function [hands] = plotIpow_selected(self, idx, fig_ax)
      hands = gobjects(1, length(idx));
      for k=1:length(idx)
        hands(k) = self.step_exps{idx(k)}.plotIpow(fig_ax);
        hold on
      end
    end
    function S = TS_dat2tex(self, varargin)
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
      p.addParameter('colormap', []);
      
      p.parse(varargin{:});
      do_color = p.Results.do_color;
      ts_vec   = p.Results.ts_vec;
      color_map = p.Results.colormap;
      
      % ------------------------------------------------------
      % ------------ Build the LaTex table -------------------
      
      if do_color
        % Build a colormap, interpolated onto the TOTAL number of settling times we have.
        if isempty(color_map)
          color_map = colormap(gca, 'jet')
        end
        x = 1:size(color_map,1);
        if isempty(ts_vec)
          ts_vec = concat_TS_cell(TS_dat_cell);
        end
        ts_vec = unique(sort(ts_vec, 'descend'));
        xq = linspace(1, (size(color_map,1)), length(ts_vec));
        
        mp_fine = zeros(length(xq), 3);
        mp_fine(:,1) = interp1(x, color_map(:,1), xq);
        mp_fine(:,2) = interp1(x, color_map(:,2), xq);
        mp_fine(:,3) = interp1(x, color_map(:,3), xq);
        
      end
      
      % -- First, we programmatically construct \tabular{ccc},
      %    since the 'ccc' depends on how many columns we need.
      %    the number of columns in the table
      [n_steps, n_exps] = size(self.TS_mat);
      
      c_fmt = repmat('c', 1, 3+n_steps);
      S = sprintf('\\begin{tabular}{%s}\n', c_fmt);
      
      % -- Form the table header:
      str_ref_cols = sprintf('&ref & delta');
      str_dat_cols = '';
      for k=1:n_exps
        str_dat_cols = sprintf(' %s & %s', str_dat_cols, self.step_exps{k}.name);
      end
      S = sprintf('%s%s%s\\\\\n\\toprule\n', S, str_ref_cols, str_dat_cols);
      %
      % -- Build up the body of the table. Outer loop is for each row.
      % Inner loop is for each experiment (columns).
      for k = 1:n_steps  %2:length(self.step_ref.step_amps)
        
        delta_ref = self.step_ref.step_diff_amps(k+1);
        str_ref_cols = sprintf('&%.2f & %.2f',...
          self.step_ref.step_amps(k+1)*self.step_ref.yscaling, delta_ref*self.step_ref.yscaling);
        str_dat_cols = '';
        % across columns, each experiment
        for j = 1:n_exps
          % settle times are stored as columns, but we have to build each row.
          ts_kj = self.TS_mat(k, j); %TS_dat_cell{j}.ts_s(k-1); % kth row, jth col.
          
          % Indexes from slowest to fastest.
          if do_color
            idx = find(ts_vec == ts_kj);
            str_dat_cols = sprintf('%s &\\cellcolor[rgb]{%.4f, %.4f, %.4f} %.2f', str_dat_cols, ...
              mp_fine(idx, 1), mp_fine(idx,2), mp_fine(idx, 3), 1000*ts_kj);
          else
            str_dat_cols = sprintf('%s & %.2f', str_dat_cols,  1000*ts_kj);
          end
        end
        s_row = sprintf('%s%s\\\\ \n', str_ref_cols, str_dat_cols);
        S = sprintf('%s%s', S, s_row);
      end
      % -- Build the footer. This is where the totals go.
      str_ref_cols = sprintf('total & -- & --');
      str_dat_cols = '';
      for j = 1:n_exps
        str_dat_cols = sprintf('%s &%.2f', str_dat_cols, ...
          1000*sum(self.TS_mat(:,j)));
      end
      
      s_row = sprintf('%s%s\\\\ \n', str_ref_cols, str_dat_cols);
      
      S = sprintf('%s\\midrule\n %s', S, s_row);
      S = sprintf('%s\\end{tabular}\n', S);
      
    end
  end % methods
  
  methods (Static)
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
    
    function write_tex_data(S, fpath)
      % write_tex_data(S, fpath)
      fid = fopen(fpath, 'w+');
      fprintf(fid, '%s', S);
      fclose(fid);
    end
    
    function ts_vec = concat_TS_cell(ts_cell)
      ts_vec = []
      for k=1:length(ts_cell)
        ts_vec = [ts_vec; ts_cell{k}.ts_s];
      end
    end
    
  end
end
