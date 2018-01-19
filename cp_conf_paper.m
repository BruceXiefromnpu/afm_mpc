% This script is only useful if you first run 
% /home/arnold/matlab/publications/afmMPC_rev1/run_results_DXP_mpc_vs_lin_exploremaxstep.m
% since it uses those workspace variables.


clc
Qx = blkdiag(Qmpc_x, zeros(Nd_x, Nd_x), 0);
Qpx = dare(sys_recyc_x.A, sys_recyc_x.B, Qx, R(1,1));

mpcp = condensedMPCprob_OA(sys_recyc_x, N_mpc, Qx, Qpx, R(1,1));
mpcp.add_U_constraint('box', slr);
a = sys_recyc_x.a(1:end-1, 1:end-1);
b = sys_recyc_x.b(1:end-1);
c = sys_recyc_x.c(1:end-1);



Nxx = SSTools.getNxNu(sys_recyc_x);
sim_struct.PLANT = ss(a,b,c,0,Ts)
% sim_struct.PLANT.InputDelay = 10;
% sim_struct.PLANT = absorbDelay(sim_struct.PLANT);
sim_struct.K_lqr = sys_recyc_x.c*0;
sim_struct.trun = Ts*800;
sim_struct.mpcProb1 = mpcp;
sim_struct.du_max = slr;

sim_struct.xss = Nxx(1:end-1)*0;
sim_struct.Nx = Nxx;
% sim_struct.x0 = Nxx
sim_struct.Ts = Ts;
sim_struct.mpc_on = 1;

[Y, U, dU] = sim_MPC_fp(sim_struct, y_ref(1)*2);

figure(300
subplot(2,1,1)
plot(Y.Time, Y.Data)

subplot(2,1,1)
plot(dU.Time, du.Data)