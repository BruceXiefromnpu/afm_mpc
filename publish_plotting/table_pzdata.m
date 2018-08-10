

addpath(fullfile(getMatPath(), 'afm_mpc_journal', 'functions'))
addpath(fullfile(getMatPath(), 'afm_mpc_journal', 'functions', 'canon'))

clear, clc
saveon = true;
plants = CanonPlants.plants_ns14();
G = plants.SYS;
G2 = plants.Gvib;


p1 = pole(G);
p2 = pole(absorbDelay(G2));

assert( all(abs(sort(p1) - sort(p2)) < 1e-13) == true)

S = ss2tex(G)

fname = fullfile(PATHS.MPCJ_root, 'latex', 'Gvib_data.tex');
fid = fopen(fname, 'w+');
fprintf(fid, '%s', S);

fclose(fid);