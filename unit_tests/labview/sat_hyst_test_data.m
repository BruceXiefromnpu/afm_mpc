
Ts = 0.01;
u = sin(2*pi*(0:200)*Ts)*10
figure, plot(u)


dp = fi(plants.hyst_sat.dp, 1, 16, 11);
wsp = fi(plants.hyst_sat.wsp, 1, 16, 11)

u_sat = PIHyst.sat_op(u, dp, wsp);
hold on
plot(u_sat)
%%

fname = 'C:\Users\arnold\Documents\matlab\afm_mpc_journal\labview\UnitTests\fpga_harnesses\Sat_Hyst_TestCase\sat_data.csv';
fid = fopen(fname, 'w+');

fprintf(fid, '%.1f\n', length(plants.hyst_sat.dp));

fprintf(fid, '%.15f, ', plants.hyst_sat.dp);
fprintf(fid, '\n');

fprintf(fid, '%.15f, ', plants.hyst_sat.wsp);
fprintf(fid, '\n');


fprintf(fid, '%.15f, ', u);
fprintf(fid, '\n');


fprintf(fid, '%.15f, ', u_sat);
fprintf(fid, '\n');

%%

rp = fi(plants.hyst_sat.rp, 1, 16, 11);
wp = fi(plants.hyst_sat.wp, 1, 16, 11)


u_hyst = PIHyst.hyst_play_op(u, rp, wp, fi(wp*0, 1, 16, 11))

figure
plot(u)
hold on
plot(u_hyst)

%%
clc

fname = 'C:\Users\arnold\Documents\matlab\afm_mpc_journal\labview\UnitTests\fpga_harnesses\Sat_Hyst_TestCase\Hyst_data.csv';
fid = fopen(fname, 'w+');

fprintf(fid, '%.1f\n', length(rp));

fprintf(fid, '%.15f, ', rp);
fprintf(fid, '\n');

fprintf(fid, '%.15f, ', wp);
fprintf(fid, '\n');


fprintf(fid, '%.15f, ', u(1:end-1));
fprintf(fid, '%.15f', u(end));
fprintf(fid, '\n');


fprintf(fid, '%.15f, ', double(u_hyst(1:end-1)));
fprintf(fid, '%.15f', double(u_hyst(end)));
fprintf(fid, '\n');
