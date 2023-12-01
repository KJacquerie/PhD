clear all
close all
clc

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

directory_name = "/Users/kathleen/Documents/PhD/2023-Project"; 
experiment_name = "Graupner2012"; 
bound_type = "SB";
fMax=50; 
new_curves=1; 
Fig_num = "Fig3"; 

Ncycles=4; 
Tcycle=40000; 
dt=0.01;
T = Ncycles*Tcycle;
Tdt = T/dt;
t = dt:dt:T;

  

w_rand= load(sprintf('%s/%s/data/%s/w.dat', directory_name, Fig_num, experiment_name));
w_STATE= load(sprintf('%s/%s/data/%s/w_BACK.dat', directory_name, Fig_num, experiment_name));

l_rand= load(sprintf('%s/%s/data/%s/late_w.dat', directory_name, Fig_num, experiment_name));
l_STATE= load(sprintf('%s/%s/data/%s/g_BACK.dat', directory_name, Fig_num, experiment_name));

wl_rand = l_rand.*w_rand;
wl_STATE = l_STATE.*w_STATE; 


    

%%
Vbin = load(sprintf('%s/%s/data/%s/Vbin.dat',directory_name, Fig_num, experiment_name));
Vspk = load(sprintf('%s/%s/data/%s/Vspk.dat',directory_name, Fig_num, experiment_name));

if(bound_type=="SB")
    V1 = load(sprintf('%s/%s/data/%s/V1.dat',directory_name, Fig_num, experiment_name));
    V2 = load(sprintf('%s/%s/data/%s/V2.dat',directory_name, Fig_num, experiment_name));
    LFP= load(sprintf('%s/%s/data/%s/LFP.dat',directory_name, Fig_num, experiment_name));
end
gion= load(sprintf('%s/%s/data/%s/gion.dat',directory_name, Fig_num, experiment_name));
gsyn= load(sprintf('%s/%s/data/%s/gsyn.dat',directory_name, Fig_num, experiment_name));

% SPB_tonic= load(sprintf('%s/%s/data/%s/SPB_tonic.dat',directory_name, Fig_num, experiment_name));
% SPB_burst= load(sprintf('%s/%s/data/%s/SPB_burst.dat',directory_name, Fig_num, experiment_name));
% PER_tonic= load(sprintf('%s/%s/data/%s/PER_tonic.dat',directory_name, Fig_num, experiment_name));
% PER_burst= load(sprintf('%s/%s/data/%s/PER_burst.dat',directory_name, Fig_num, experiment_name));
% IBF_tonic= load(sprintf('%s/%s/data/%s/IBF_tonic.dat',directory_name, Fig_num, experiment_name));
% IBF_burst= load(sprintf('%s/%s/data/%s/IBF_burst.dat',directory_name, Fig_num, experiment_name));
% DC_tonic= load(sprintf('%s/%s/data/%s/DC_tonic.dat',directory_name, Fig_num, experiment_name));
% DC_burst= load(sprintf('%s/%s/data/%s/DC_burst.dat',directory_name, Fig_num, experiment_name));
% freq_tonic= load(sprintf('%s/%s/data/%s/freq_tonic.dat',directory_name, Fig_num, experiment_name));
% freq_burst= load(sprintf('%s/%s/data/%s/freq_burst.dat',directory_name, Fig_num, experiment_name));
%  

%%


color_tonic = [253 236 214]./255; 
color_burst = [157 195 236]./255; 
color_gray = [ 0.3 0.3 0.3]; 
color_p = [ 0 0 0;0.2 0.2 0.2;0.5 0.5 0.5; 0.8 0.8 0.8];
color_p = [ 0.8 0.8 0.8;0.5 0.5 0.5; 0.2 0.2 0.2;0 0 0];


% couleur rainbow
%color_p = [200 52 93; 255 193 0;112 173 71; 237 125 4991 155 213 ]./255;


ptx_long = 5; 
ptx_short = ptx_long/2; 
pty = 2; 

ptx_hist = 3; 
pty_hist= 3;

color_hist = [0.7 0.7 0.7]; 


%% w
points = [1 19 32 47 80] ;
figure
hold on
for idx_t=1:1:4
    plot([idx_t*20e3 idx_t*20e3], [-1 1], 'linewidth', 0.25, 'color', [0.3 0.3 0.3]) 
end
%plot(w_rand', 'color', [0.75 0.75 0.75], 'linewidth', 0.25)
for idx_p=1:1:size(w_rand,1)
    plot( w_rand(idx_p,:)','color', color_p(idx_p,:),'linewidth', 1.5)
end
    

%xlim([T1 T2])
mm = min(min(w_rand)); 
MM = max(max(w_rand)); 
ylim([mm MM]); 

xticks([0 20000 40000 60000 80000 100000 120000 140000 160000])
xticklabels({'','','','','','','','',''})
yticks([mm MM])
yticklabels({'', ''})


set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 ptx_long pty]);
print(sprintf('%s/%s/fig/%s_w',directory_name, Fig_num, experiment_name), '-depsc', '-painters')


mm_MM = [mm MM]; 
csvwrite(sprintf('%s/%s/fig/mm_MM',directory_name, Fig_num), mm_MM)

%%

figure
hold on
for idx_t=1:1:4
    plot([idx_t*20e3 idx_t*20e3], [-1 1], 'linewidth', 0.25, 'color', [0.3 0.3 0.3]) 
end
%plot(l_rand', 'color', [0.75 0.75 0.75], 'linewidth', 0.25)
for idx_p=1:1:size(w_rand,1)
    plot( l_rand(idx_p,:)','color', color_p(idx_p,:),'linewidth', 1.5)
end 

%xlim([T1 T2])
mm = min(min(l_rand)); 
MM = max(max(l_rand)); 
if(mm~=MM)
    ylim([mm MM]); 
    yticks([mm MM])
    yticklabels({'', ''})
else
    ylim([mm*0.9 1.1*mm])
    yticks([mm ])
    yticklabels({''})   
end
xticks([0 20000 40000 60000 80000 100000 120000 140000 160000])
xticklabels({'','','','','','','','',''})



set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 ptx_long pty]);
print(sprintf('%s/%s/fig/%s_l',directory_name, Fig_num,experiment_name), '-depsc', '-painters')


mm_MM = [mm MM]; 
csvwrite(sprintf('%s/%s/fig/l_mm_MM',directory_name, Fig_num), mm_MM)


%%
figure
hold on
for idx_t=1:1:4
    plot([idx_t*20e3 idx_t*20e3], [-1 1], 'linewidth', 0.25, 'color', [0.3 0.3 0.3]) 
end
%plot(wl_rand', 'color', [0.75 0.75 0.75], 'linewidth', 0.25)
for idx_p=1:1:size(w_rand,1)
    plot( wl_rand(idx_p,:)','color', color_p(idx_p,:),'linewidth', 1.5)
end
    

%xlim([T1 T2])
mm = min(min(wl_rand)); 
MM = max(max(wl_rand)); 
ylim([mm MM]); 
xticks([0 20000 40000 60000 80000 100000 120000 140000 160000])
xticklabels({'','','','','','','','',''})
yticks([mm MM])
yticklabels({'', ''})


set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 ptx_long pty]);
print(sprintf('%s/%s/fig/%s_wl',directory_name, Fig_num,experiment_name), '-depsc', '-painters')


mm_MM = [mm MM]; 
csvwrite(sprintf('%s/%s/fig/wl_mm_MM',directory_name, Fig_num), mm_MM)



%% Raster
T1 = Tcycle/2-1500; 
T2 = Tcycle/2+2500;

figure
count=1;
idx_t = find(Vbin(1,:)==1); 
hold on
for idx=2:1:size(Vbin,1) % start at 2 to remove the inhibitor
    idx_t = find(Vbin(idx,:)==1); 
    
    plot(idx_t, Vbin(idx,idx_t)+count,'o','color', 'k', 'MarkerSize', 0.25)
    count=count-0.01; 
end

xlim([T1 T2])
ylim([1 2])
xticks([0 1000])
xticklabels({'',''})
yticks([0 1])
yticklabels({'',''})

axis off
box off


set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 5 5]);
print(sprintf('%s/%s/fig/Raster_%s', directory_name, Fig_num, experiment_name), '-depsc', '-painters')

%% LFP

interv = T1/dt:T2/dt;

yy=smooth(LFP,1500);

figure
hold on
plot([0 interv(end)], [0.1 0.1], '-', 'color', [0.7 0.7 0.7], 'linewidth', 0.25) 
plot([0 interv(end)], [0.2 0.2], '-', 'color', [0.7 0.7 0.7], 'linewidth', 0.25) 
plot([0 interv(end)], [0.3 0.3], '-', 'color', [0.7 0.7 0.7], 'linewidth', 0.25) 
plot([0 interv(end)], [0.4 0.4], '-', 'color', [0.7 0.7 0.7], 'linewidth', 0.25) 

plot(t(interv)*1e-3, yy(interv), 'color', color_gray, 'linewidth', 0.75) 
ylim([-0.001 0.03])
xlim([(T1-100)/1000 (T2+100)/1000])
yticks([0 0.1])
yticklabels({'',''})
xticks(0)
xticklabels({''})

box off
%axis off
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 ptx_long 1.5]);
print(sprintf('%s/%s/fig/LFP_%s', directory_name, Fig_num, experiment_name), '-depsc')

