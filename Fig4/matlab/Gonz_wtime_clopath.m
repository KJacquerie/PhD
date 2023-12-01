clear all 
close all 
clc

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

directory_name = "/Users/kathleen/Documents/PhD/2023-Fast/Gonzalez";
%%
pt = 11; 
ptx = 5; 
pty=2.5; 


SD = 0;
idx_ntk=3;

if(SD==1)
     experiment_mat = ["AMPA_SD", "Reset_SD"]; 
else
     experiment_mat = ["AMPA"]% ,"Reset"]; 
end

switch idx_ntk
    case 1
        ncellsC = 11;
    case 2
        ncellsC = 21;
    case 3
        ncellsC = 101; 
end



%%
color_blue = [106 153 208]./255; 
color_learning = [222 235 247]./255; 
color_gray = [ 0.3 0.3 0.3]; 
color_pink = [250/255 244/255 247/255];

color_corr = [0.25 0.25 0.25]; %[63/255 92/255 206/255]; 
color_uncorr = [0.75 0.75 0.75];%[130/255 187/255 255/255]; 
color_green = [112 173 71]./255; 
color_reset = [228 162 180]./255; %[ 157 12 58]./255; 

%%
for idx_expm=1:1:length(experiment_mat)

    experiment_name = experiment_mat(idx_expm); 
w = load(sprintf('%s/data/%s_w%d.dat', directory_name, experiment_name, idx_ntk));
g = load(sprintf('%s/data/%s_g%d.dat', directory_name,experiment_name, idx_ntk));

%%


connections = [1, 19, 30, 69, 81];

%%
figure
if(SD==1)
    bgd_SD
else
    bgd_seq; 
end

    hold on
    plot(w', 'color', color_uncorr, 'linewidth', 1)
    ylim([0 1])

    for i =1:1:length(connections)
        hold on
        plot(w(connections(i),:), 'color', color_corr, 'linewidth', 1)
        ylim([0 1])
    end
    
xlim([0 150000])
%axis off
ylim([0 1])
box off

set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 ptx pty]);
print(sprintf('%s/fig_clopath/%s_w%d',directory_name, experiment_name, idx_ntk),'-depsc', '-painters')


%%
figure
if(SD==1)
    bgd_SD
else
    bgd_seq; 
end

    hold on
    plot(g', 'color', color_uncorr, 'linewidth', 1)
    ylim([0 1])

    for i =1:1:length(connections)
        hold on
        plot(g(connections(i),:), 'color', color_corr, 'linewidth', 1)
        ylim([0 1])
    end
    

box off
ylim([0 0.0252])
%ylim([0 1.1*max(max(g))])
%axis off
xlim([0 150000])
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 ptx pty]);
print(sprintf('%s/fig_clopath/%s_g%d', directory_name,experiment_name, idx_ntk),'-depsc', '-painters')


%%

%subplot(3,1,3)
figure    
if(SD==1)
    bgd_SD
else
    bgd_seq; 
end
    hold on
    plot(w'.*g', 'color', color_uncorr, 'linewidth', 1)
    ylim([0 1])

    for i =1:1:length(connections)
        hold on
        plot(w(connections(i),:).*g(connections(i),:), 'color', color_corr, 'linewidth', 1)
        ylim([0 1])
    end   
    
xlim([0 150000])
%ylim([0 max(max(w.*g))])
ylim([0 0.0162])
%axis off 
box off
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 ptx pty]);
print(sprintf('%s/fig_clopath/%s_wg%d', directory_name,experiment_name, idx_ntk),'-depsc', '-painters')


%%
for i=1:1:length(connections)

    wsave(i+(idx_ntk-1)*10,:) = w(connections(i),:); 
end

for i=1:1:length(connections)
    gsave(i+(idx_ntk-1)*10,:) = g(connections(i),:);
end


    T_wake  = [14999, 44999, 74999, 104999, 134999];
    T_sleep = [29999, 59999, 89999, 119999, 149999];
%%
    gtot = g.*w;
    SNR_wake = zeros(1,length(T_wake)+1);
    SNR_sleep = zeros(1,length(T_sleep)+1);
    SNR_test = 0; 

    
    
for j=1:1:length(T_wake)
   % for i=1:1:length(connections)
        w_wake(:,j) = gtot(:,T_wake(j));
        w_sleep(:,j) = gtot(:,T_sleep(j));
        w_test(:,j) = gtot(:,end);
    %end

    SNR_wake(j)  = max(w_wake(:,j))/mean(w_wake(:,j)); 
    SNR_sleep(j) = max(w_sleep(:,j))/mean(w_sleep(:,j)); 
    SNR_test = max(w_test(:,j))/mean(w_test(:,j)); 
end


SNR = [SNR_wake; SNR_sleep; 0 0 0 0 0 SNR_test]; 

%%

figure
%bar(SNR')
hold on
bar([1 2 3 4 5], SNR_wake(1:5),'FaceColor',color_blue,'EdgeColor',color_blue,'LineWidth',0.5, 'BarWidth', 0.25)
if(SD==1)
    bar([1.5 2.5 3.5 4.5 5.5], SNR_sleep(1:5),'FaceColor','white','EdgeColor',color_gray,'LineWidth',0.5, 'BarWidth', 0.25)
else
    bar([1.5 2.5 3.5 4.5 5.5], SNR_sleep(1:5),'FaceColor',color_reset,'EdgeColor',color_reset,'LineWidth',0.5, 'BarWidth', 0.25)

end


xlim([0.75 5.75])
ylim([0 21])

%axis off
box off

box off
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 ptx pty*2]);
print(sprintf('%s/fig_clopath/%s_SNR%d', directory_name,experiment_name, idx_ntk),'-depsc', '-painters')


end
    
