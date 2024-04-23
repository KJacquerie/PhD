clear all
close all
clc

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

directory_name = "/Users/kathleen/Documents/PhD/2023-Project"; 
experiment_name = "Graupner2012_TUNE"; 
bound_type = "SB";
fMax=50; 
new_curves=1; 
Fig_num = "Fig2"; 

Ncycles=1; 
Tcycle=80000; 
Nstates=4; 
dt=0.01;
T = Ncycles*Tcycle;
Tdt = T/dt;
t = dt:dt:T;

%tmp =  load(sprintf('%s/%s/data/%s/idx_w.mat', directory_name, Fig_num, experiment_name));
%idx_w = tmp.idx_w;
idx_w = load('idx_wl.dat');


  


  w_rand= load(sprintf('%s/%s/data/%s/w.dat', directory_name, Fig_num, experiment_name));
  %wtemp = w_rand(idx_wl, :);
  %w_rand =wtemp; 
  w_STATE= load(sprintf('%s/%s/data/%s/w_BACK.dat', directory_name, Fig_num, experiment_name));
  
 %%
Vbin = load(sprintf('%s/%s/data/%s/Vbin.dat',directory_name, Fig_num, experiment_name));
Vspk = load(sprintf('%s/%s/data/%s/Vspk.dat',directory_name, Fig_num, experiment_name));

if(bound_type=="SB")
    %V1 = load(sprintf('%s/%s/data/%s/V1.dat',directory_name, Fig_num, experiment_name));
    %V2 = load(sprintf('%s/%s/data/%s/V2.dat',directory_name, Fig_num, experiment_name));
    LFP= load(sprintf('%s/%s/data/%s/LFP.dat',directory_name, Fig_num, experiment_name));
end
gion= load(sprintf('%s/%s/data/%s/gion.dat',directory_name, Fig_num, experiment_name));
gsyn= load(sprintf('%s/%s/data/%s/gsyn.dat',directory_name, Fig_num, experiment_name));

SPB_tonic= load(sprintf('%s/%s/data/%s/SPB_tonic.dat',directory_name, Fig_num, experiment_name));
SPB_burst= load(sprintf('%s/%s/data/%s/SPB_burst.dat',directory_name, Fig_num, experiment_name));
PER_tonic= load(sprintf('%s/%s/data/%s/PER_tonic.dat',directory_name, Fig_num, experiment_name));
PER_burst= load(sprintf('%s/%s/data/%s/PER_burst.dat',directory_name, Fig_num, experiment_name));
IBF_tonic= load(sprintf('%s/%s/data/%s/IBF_tonic.dat',directory_name, Fig_num, experiment_name));
IBF_burst= load(sprintf('%s/%s/data/%s/IBF_burst.dat',directory_name, Fig_num, experiment_name));
DC_tonic= load(sprintf('%s/%s/data/%s/DC_tonic.dat',directory_name, Fig_num, experiment_name));
DC_burst= load(sprintf('%s/%s/data/%s/DC_burst.dat',directory_name, Fig_num, experiment_name));
freq_tonic= load(sprintf('%s/%s/data/%s/freq_tonic.dat',directory_name, Fig_num, experiment_name));
freq_burst= load(sprintf('%s/%s/data/%s/freq_burst.dat',directory_name, Fig_num, experiment_name));

%%


color_tonic = [253 236 214]./255; 
color_burst = [157 195 236]./255; 
color_gray = [ 0.3 0.3 0.3]; 

% couleur rainbow
color_p = [200 52 93; 
            255 193 0;
            112 173 71;
            237 125 49
            91 155 213; 
            ]./255;


ptx_long = 15; 
ptx_short = ptx_long/2; 
pty = 4; 

ptx_hist = 3; 
pty_hist= 3;

color_hist = [0.7 0.7 0.7]; 


%% w
points = [3 19 32 47 100] ;
figure
hold on
%bgd_seq(Ncycles, Tcycle, color_tonic,color_burst)
%if(new_curves==1)
%    points = [10 100 1000 2000 2500] ;
%    idx_w =randi([1 size(w,1)],100,1);
%    plot(w(idx_w,:)', 'color', [0.75 0.75 0.75], 'linewidth', 0.25)
%    for idx_p=1:1:length(points)
%        plot( w(points(idx_p),:)','color', color_p(idx_p, :),'linewidth', 2)
%    end
%else

for idx_t=1:1:4
    plot([idx_t*20e3 idx_t*20e3], [0 1], 'linewidth', 0.25, 'color', [0.3 0.3 0.3]) 
end
plot(w_rand', 'color', [0.75 0.75 0.75], 'linewidth', 0.25)
for idx_p=1:1:length(points)
    plot( w_rand(points(idx_p),:)','color', color_p(idx_p, :),'linewidth', 2)
end
    
mm = min(min(w_STATE)); 
MM = max(max(w_STATE)); 
ylim([mm MM]); 
xticks([0 20000 40000 60000 80000 100000 120000 140000 160000])
xticklabels({'','','','','','','','',''})
yticks([mm MM])
yticklabels({'', ''})


set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 ptx_long pty*0.6]);
print(sprintf('%s/%s/fig/%s/w',directory_name,  Fig_num, experiment_name), '-depsc', '-painters')


mm_MM = [mm MM]; 
csvwrite(sprintf('%s/%s/fig/%s/w_mm_MM',directory_name,  Fig_num, experiment_name), mm_MM)

%% Raster
Tstate=20000; 
for idx_state=1:2:8
    T1 = idx_state*Tstate/2-1000; 
    T2 = idx_state*Tstate/2+1000;

    figure
    count=1;
    idx_t = find(Vbin(1,:)==1); 
    hold on
    for idx=2:1:size(Vbin,1) % start at 2 to remove the inhibitor
        idx_t = find(Vbin(idx,:)==1); 

        plot(idx_t, Vbin(idx,idx_t)+count,'o','color', 'k', 'MarkerSize', 0.15)
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
    set(gcf, 'PaperPosition', [0 0 3 2.5]);
   print(sprintf('%s/%s/fig/%s/Raster%d',directory_name,  Fig_num, experiment_name, idx_state), '-depsc', '-painters')

    
    % _____ LFP

    T1 = idx_state*Tstate/2-1000; 
    T2 = idx_state*Tstate/2+1000;


interv = T1/dt:T2/dt;

yy=smooth(LFP,1500);

figure
plot(t(interv)*1e-3, yy(interv), 'color', color_gray, 'linewidth', 0.5) 
ylim([-0.05 0.4])
xlim([(T1-100)/1000 (T2+100)/1000])
yticks([0 0.1])
yticklabels({'',''})
xticks(0)
xticklabels({''})

box off
%axis off
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 3.5 1.5]);
print(sprintf('%s/%s/fig/%s/LFP%d',directory_name,  Fig_num, experiment_name, idx_state), '-depsc', '-painters')
end



%% HIST

figure
for idx=1:1:size(SPB_burst,2)
    histogram(SPB_burst(:,idx),'BinEdges', 0.5:1:5.5, 'FaceColor', color_hist)
    xticks([1 2 3 4 5])
    xticklabels({'','','','',''})
    yticks([0  100])
    yticklabels({'',''})
    ylim([0 100])
    xlim([1 6])
    box off
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 ptx_hist pty_hist]);
    print(sprintf('%s/%s/fig/%s/SPB%d', directory_name, Fig_num,experiment_name, idx), '-depsc', '-painters')
end
%%

figure
for idx=1:1:size(IBF_burst,2)
    histogram(IBF_burst(:,idx),'BinEdges', 0:10:100,  'FaceColor', color_hist)
    xticks([0 100])
    xticklabels({'',''})
    yticks([0 100])
    yticklabels({'',''})
    ylim([0 100])
    xlim([0 100])
    box off
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 ptx_hist pty_hist]);
    print(sprintf('%s/%s/fig/%s/IBF%d',directory_name, Fig_num,experiment_name,idx), '-depsc', '-painters')
end

%%

for idx=1:1:size(PER_burst,2)
    figure(idx)
    histogram(1000./PER_burst(:,idx),'BinEdges', 2:0.05:8, 'FaceColor', color_hist)
    xticks([0])
    xticklabels({''})
    yticks([0  100])
    yticklabels({'',''})
    ylim([0 100])
    xlim([2 8])
    box off
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 ptx_hist pty_hist]);
    print(sprintf('%s/%s/fig/%s/PER%d',directory_name, Fig_num,experiment_name,idx), '-depsc', '-painters')
end
%%

figure
for idx=1:1:size(freq_tonic,2)
    histogram(freq_tonic(:,idx),'BinEdges', 0:10:70, 'FaceColor', color_hist)
    xticks([0 70])
    xticklabels({'',''})
    yticks([0  100])
    yticklabels({'',''})
    ylim([0 100])
    xlim([0 70])
    box off
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 ptx_hist pty_hist]);
    print(sprintf('%s/%s/fig/%s/freq%d',directory_name, Fig_num,experiment_name,idx), '-depsc', '-painters')
end

%%

edges = [0:0.025:1];
for idx=1:1:size(w_STATE,2)
    figure
    histogram(w_STATE(:,idx)','BinEdges', edges,'FaceColor', color_hist)
    xticks([0 1])
    xticklabels({'',''})
    yticks([0 1000])
    yticklabels({'',''})
    ylim([-100 1250])
    xlim([0 1])
    box off
    %axis off
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 ptx_hist*1.5 pty_hist*0.8]);
    print(sprintf('%s/%s/fig/%s/whist%d',directory_name, Fig_num,experiment_name, idx), '-depsc', '-painters')
end

%% 

edges = [0:0.025:1];
color_H = [0.7 0.7 0.7;  0 0 0.2; 0 0 0.5]
figure
for idx=2:1:4%size(w_STATE,2)
    hold on
    histogram(w_STATE(:,idx)','BinEdges', edges,'FaceColor', color_H(idx-1,:))
    xticks([0 1])
    xticklabels({'',''})
    yticks([0 1000 ])
    yticklabels({'',''})
    ylim([-100 1250])
    xlim([0 1])
    box off
    %axis off
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 ptx_hist*1.5 pty_hist*0.8]);
    print(sprintf('%s/%s/fig/%s/whistSUP',directory_name, Fig_num,experiment_name), '-depsc', '-painters')
end


