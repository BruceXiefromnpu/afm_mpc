% saveControlData(AllMatrix Ki, Nbar, Ns, Nd, r, y_uKx, controlDataPath, refTrajPath)
%
% Save controller data and reference trajectories to .csv files for reading
% by labview.
%
%   Inputs:
%       AllMatrix: [B; L; K'; C'; A(1,:)....A(Ns,:)]. Create with
%          packMatrix(sys, L, K)
%       Ki: integral gain, scalar
%       Nbar: scalar, ff gain
%       Ns: scaler, # states
%       Nd: scalar, # delay states
%       r: scaler, set-point (ref)
%
%       The previous inputs are saved into a tall, 2 column matrx:
%           [AllMatrix, [Ki; Nbar; Ns;Nd;r] ], to controlDataPath
%       
%       y_uKx: vector of interleaved data. Genereate with 
%           pack_uKx_y(uKx, y, t)
%       controlDataPath: string, path to save controlData vector.
%       refTrajPath: string, path to save y_uKx vector.

function saveControlData(y_uKx, refTrajPath)

dlmwrite(refTrajPath, y_uKx, 'delimiter', ',', 'precision', 12)         


end