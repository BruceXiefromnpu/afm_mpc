clear 

load hyst_sat_data.mat

whos

fid = fopen('hyst_sat_data.csv', 'w+');

% fprintf(fid, 'wp, %f\n', hyst_sat_data.wp')
fprintf(fid, 'Nyst, %d\n', length(hyst_sat_data.wsp));
fprintf(fid, 'Nsat, %d\n', length(hyst_sat_data.wp));
fprintf(fid, 'Ntraj, %d\n', length(hyst_sat_data.u_gdinv.Data));

write_labeled_row(fid, 'wp', hyst_sat_data.wp);
write_labeled_row(fid, 'rp', hyst_sat_data.rp);
write_labeled_row(fid, 'wsp', hyst_sat_data.wsp);
write_labeled_row(fid, 'dp', hyst_sat_data.dp);

write_labeled_row(fid, 'u_gdinv', hyst_sat_data.u_gdinv.Data);
write_labeled_row(fid, 'u_sinv', hyst_sat_data.u_sinv.Data);
write_labeled_row(fid, 'u_hinv', hyst_sat_data.u_hinv.Data);




fclose(fid);
%%
wp = hyst_sat_data.wp;
rp = hyst_sat_data.rp;
PIHyst.hyst_play_op(hyst_sat_data.u_sinv.Data(1:2), rp, wp, wp*0)


%%

function write_labeled_row(fid, label, rowdata)
  
  fprintf(fid, '%s,', label);  
  fprintf(fid, '%f,', rowdata(1:end-1));
  fprintf(fid, '%f', rowdata(end));
  
  fprintf(fid, '\n');
  
end
