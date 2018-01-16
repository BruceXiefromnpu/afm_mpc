function [Y, U, dU] = sim_MPC_fp(sim_struct, ref_f)
% ------------------------------------------------------------------- %
% Pull out all the data stored in sim_struct to expose it to
% simulink. There must be a better way...

    options = simset('SrcWorkspace','current');
    % Expose the sime struct to simulink.
    K_lqr = sim_struct.K_lqr;
    PLANT = sim_struct.PLANT;
    trun = sim_struct.trun;
    mpcProb1 = sim_struct.mpcProb1;
    du_max = sim_struct.du_max;
    mpc_on = sim_struct.mpc_on;
    xss = sim_struct.xss;
    Nx = sim_struct.Nx;
    
    x0 = sim_struct.xss*0;
    uss_0 = 0;
    Ts = PLANT.Ts;
    sim('MPC_fp', [], options)
    
    Y = y_mpcDist;
    U = u_mpcDist;
    dU = du_mpcDist;
end
