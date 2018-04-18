

MF_unbolted = load(['/media/labserver/mpc-journal/sysID/x-' ...
                    'axis_sines_info_out_3-30-2018-01.mat']);
MF_bolted = load(['/media/labserver/mpc-journal/sysID/x-' ...
                    'axis_sines_info_out_3-30-2018-03.mat']);

G_unbolted = MF_unbolted.modelFit.frf.G_frf(:);
freqs_unbolted = MF_unbolted.modelFit.frf.freq_s;


G_bolted = MF_bolted.modelFit.frf.G_frf(:);
freqs_bolted = MF_bolted.modelFit.frf.freq_s;

F1 = figure(1);

h2 = frfBode(G_unbolted, freqs_unbolted, F1, '-b', 'Hz');

h1 = frfBode(G_bolted, freqs_bolted, F1, '--r', 'Hz');


h1.DisplayName = 'Stage Bolted to AFM base';
h2.DisplayName = 'Stage free on table';

subplot(2,1,1)
leg1 = legend([h1, h2]);
set(leg1, 'location', 'SouthWest')
F1.PaperPosition = [1.3333    2.2135    5.8333    6.5729];

if 0
    F1.PaperPosition = [1.3333    2.2135    5.8333    6.5729];
    saveas(F1, '/home/arnold/matlab/afm_mpc_journal/figures/pow-amp/unbolted_frf.svg')
end
