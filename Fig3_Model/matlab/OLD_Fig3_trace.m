clear all
close all
clc

%%%% ALLER DANS OLD TRACE %%%%%
% car il n'y a plus l'enregistrement de w.dat
%
%
%
%
%
%
%
%
%
%
%
%



set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

directory_name = "/Users/kathleen/Documents/PhD/2023-reset"; 
experiment_name = "Graupner2012"; 
region = "Cortex"; 
bound_type = "SB";

Fig_num = "Fig2"; 


Vbin = load(sprintf('%s/%s/data/%s/%s/%s/Vbin.dat', directory_name,Fig_num, experiment_name, region, bound_type));
Vspk = load(sprintf('%s/%s/data/%s/%s/%s/Vspk.dat', directory_name,Fig_num, experiment_name, region, bound_type));
if(plot_w==1)
    tmp =  load(sprintf('%s/Fig1/data/%s/%s/idx_w.mat', directory_name, experiment_name, bound_type));
    idx_w = tmp.idx_w;
    w = load(sprintf('%s/%s/data/%s/%s/%s/w.dat',  directory_name,Fig_num, experiment_name, region, bound_type));
    w_rand = w(idx_w,:); 
    save(sprintf('%s/Fig2/data/%s/%s/w_rand.mat', directory_name, experiment_name, bound_type), 'w_rand');
  
end
w_BACK = load(sprintf('%s/%s/data/%s/%s/%s/w.dat',  directory_name,Fig_num, experiment_name, region, bound_type));
tmp =  load(sprintf('%s/Fig1/data/%s/%s/idx_w.mat', directory_name, experiment_name, bound_type));
idx_w = tmp.idx_w;
 
gion= load(sprintf('%s/%s/data/gion.dat', directory_name,Fig_num)); 
gsyn= load(sprintf('%s/%s/data/gsyn.dat', directory_name,Fig_num));
%%
% SPB_tonic= load(sprintf('%s/%s/data/SPB_tonic.dat', directory_name,Fig_num));
% SPB_burst= load(sprintf('%s/%s/data/SPB_burst.dat', directory_name,Fig_num));
% PER_tonic= load(sprintf('%s/%s/data/PER_tonic.dat', directory_name,Fig_num));
% PER_burst= load(sprintf('%s/%s/data/PER_burst.dat', directory_name,Fig_num));
% IBF_tonic= load(sprintf('%s/%s/data/IBF_tonic.dat', directory_name,Fig_num));
% IBF_burst= load(sprintf('%s/%s/data/IBF_burst.dat', directory_name,Fig_num));
% DC_tonic= load(sprintf('%s/%s/data/DC_tonic.dat', directory_name,Fig_num));
% DC_burst= load(sprintf('%s/%s/data/DC_burst.dat', directory_name,Fig_num));
% freq_tonic= load(sprintf('%s/%s/data/freq_tonic.dat', directory_name,Fig_num));
% freq_burst= load(sprintf('%s/%s/data/freq_burst.dat', directory_name,Fig_num));

%%



Ncycles=1; 
Tcycle=40000; 
dt=0.01;
T = Ncycles*Tcycle;
Tdt = T/dt;
t = dt:dt:T;


color_tonic = [253 236 214]./255; 
color_burst = [157 195 236]./255; 
color_gray = [ 0.3 0.3 0.3]; 


ptx_long = 15; 
ptx_short = ptx_long/2; 
pty = 3; 
ptx=3;

%% w

figure
hold on
%bgd_seq(Ncycles, Tcycle, color_tonic,color_burst)
load('idx_w.mat');
load('idx_choice.mat');
plot(w(idx_w,:)', 'color', [0.75 0.75 0.75], 'linewidth', 0.25)
if(bound_type=="SB")
    plot(w(idx_w(20),:), 'color', 'k', 'linewidth', 1.5)
    ylim([0 1])
else
    plot(w(idx_choice,:), 'color', 'k', 'linewidth', 1.5)
end

xticks([0 20000 40000 ])
xticklabels({'','',''})
yticks([0 0.5 1])
yticklabels({'','', ''})


set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 ptx pty]);
print(sprintf('/Users/kathleen/Documents/PhD/2023-reset/%s/fig/trace/%s_%s_%s_w',Fig_num, experiment_name, region, bound_type), '-depsc', '-painters')

%%