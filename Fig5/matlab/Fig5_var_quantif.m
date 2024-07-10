clear all
close all
clc

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

directory_name = "/Users/kathleen/Documents/PhD/2023-Project"; 
experiment_name = "Graupner2012_VAR"; 
experiment_type = "normal"; 
bound_type = "SB";
fMax=50; 
new_curves=0; 
Fig_num = "Fig5"; 

tau_mat = [10 20 30 40 50 60 70 80 90 100 200 300 400 500 600 700 800 900 1000 2000] % 100:100:1000;
mat = [0.0001 0.00025 0.0005 0.00075 0.001 0.0025 0.005 0.0075 0.01];
NB_l = length(mat);
SNR_count=zeros(length(tau_mat), NB_l);

%tau_val= 100;

for idx_tau = 1:1:size(tau_mat,2)
    tau_val = tau_mat(idx_tau);
for idx_lINIT=1:1:NB_l
    %idx_lINIT
    Ncycles=5; 
    Tcycle=30000; 
    dt=0.01;
    T = Ncycles*Tcycle;
    Tdt = T/dt;
    t = dt:dt:T;



    %w_rand= load(sprintf('%s/%s/data/%s/tau=%d/%d/w.dat', directory_name, Fig_num, experiment_name, tau_val, idx_lINIT));
    w_STATE= load(sprintf('%s/%s/data/%s/%s/tau=%d/%d/w_BACK.dat', directory_name, Fig_num, experiment_name, experiment_type, tau_val, idx_lINIT));

    %l_rand= load(sprintf('%s/%s/data/%s/tau=%d/%d/late_w.dat', directory_name, Fig_num, experiment_name, tau_val, idx_lINIT));
    l_STATE= load(sprintf('%s/%s/data/%s/%s/tau=%d/%d/g_BACK.dat', directory_name, Fig_num, experiment_name, experiment_type, tau_val, idx_lINIT));

    %wl_rand = l_rand.*w_rand;
    wl_STATE = l_STATE.*w_STATE; 





    %%


    ptx_long = 5; 
    ptx_short = ptx_long/2; 
    pty = 2; 

    ptx_hist = 3; 
    pty_hist= 3;

    color_hist = [0.7 0.7 0.7]; 



    %%

    color_tonic = [31 78 121]./255; %
    color_burst = [189 215 238]./255; 
    color_blue= [106 153 208]./255; 
    color_learning = [222 235 247]./255; 
    color_gray = [ 0.3 0.3 0.3]; 
    color_pink = [250/255 244/255 247/255];

    color_corr = [0.25 0.25 0.25]; %[63/255 92/255 206/255]; 
    color_uncorr = [0.75 0.75 0.75];%[130/255 187/255 255/255]; 
    color_green = [112 173 71]./255; 
    color_reset = [228 162 180]./255; %[ 157 12 58]./255; 

    

    %%
    count=1; 
    for idx_state=1:2:Ncycles*2-1

        wl_tonic(:,count) = wl_STATE(:,idx_state);
        wl_burst(:,count) = wl_STATE(:,idx_state+1);

        SNR_tonic(count)  = max(wl_tonic(:,count))/mean(wl_tonic(:,count)); 
        SNR_burst(count) = max(wl_burst(:,count))/mean(wl_burst(:,count)); 
        
        if(SNR_tonic(count)>SNR_burst(count))
            SNR_count(idx_tau, idx_lINIT) = SNR_count(idx_tau, idx_lINIT) -1; % downselection
            SNR_param{idx_tau,idx_lINIT}(count)=-1;
        else
            if(SNR_tonic(count)<SNR_burst(count))
                SNR_count(idx_tau, idx_lINIT) = SNR_count(idx_tau, idx_lINIT) +1; % up selection
                SNR_param{idx_tau,idx_lINIT}(count)=+1;
            else
                print('here')
            end
        end
        %SNR_tonic(count)  = mean(wl_tonic(1:5,count))/mean(wl_tonic(6:100,count)); 
        %SNR_burst(count) = mean(wl_burst(1:5,count))/mean(wl_burst(6:100,count)); 
        count=count+1; 
        
        SNR_tonic_MAXX{idx_tau, idx_lINIT} = max(SNR_tonic); 
        SNR_burst_MAXX{idx_tau, idx_lINIT} = max(SNR_burst); 
    end


    SNR = [SNR_tonic; SNR_burst];

    figure
    %bar(SNR')
    %subplot(4,1,4)
    hold on
    bar([1:2:2*size(SNR_tonic,2)], SNR_tonic(1:end),'FaceColor',color_tonic,'EdgeColor',color_tonic,'LineWidth',0.5, 'BarWidth', 0.25)
    bar([2:2:2*size(SNR_tonic,2)], SNR_burst(1:end),'FaceColor',color_burst,'EdgeColor',color_burst,'LineWidth',0.5, 'BarWidth', 0.25)
    xlim([0 11])
   ylim([0 22.5])
   % ylim([0 2.5])
    box off
    %axis off
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 7 4]);
    print(sprintf('%s/%s/fig/%s/%s/SNR_%d_%d', directory_name, Fig_num,experiment_name,experiment_type,  tau_val, idx_lINIT), '-dsvg')


    close all
    
    %upselection 2
    %downselection -2
    % no change 0
    % from down to up 1
    % from up to down -1
    if(~isempty(sort([strfind(SNR_param{idx_tau,idx_lINIT}>0, [0 1]) strfind(SNR_param{idx_tau,idx_lINIT}>=0, [1 0])])))
        if(SNR_param{idx_tau,idx_lINIT}(1)>0) % from up-selection to down-selection
            SNR_val(idx_tau, idx_lINIT) = -1;
        else
            SNR_val(idx_tau, idx_lINIT) = 1;
        end
    else
        if(SNR_param{idx_tau,idx_lINIT}(1)>0)
            SNR_val(idx_tau, idx_lINIT) = 2;
        else
            SNR_val(idx_tau, idx_lINIT) = -2;
        end
    end
end
end

%% 
figure
imagesc(SNR_val)
colorbar
box off
%axis off
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 9 8]);
print(sprintf('%s/%s/fig/%s/%s/_map_val', directory_name, Fig_num,experiment_name,experiment_type), '-dsvg')


figure
imagesc(SNR_count)
colorbar
box off
%axis off
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 9 8]);
print(sprintf('%s/%s/fig/%s/%s/_map_count', directory_name, Fig_num,experiment_name,experiment_type), '-dsvg')




%%


idx_row= [1:19];  

figure
imagesc(SNR_val(idx_row,:))



hexMap = {'1E4A76', '5FA3BF', 'E17438', 'F3B15D',  '70AD47'};
hexMap = {'30708E', '7DC7D2', 'EDECED', 'E58331', 'C43424'} ;
%hexMap = {'C0C0C0', '808080', '404040', '000000', 'FF99CC', '9999FF', '3333FF', '000099', '3399FF', '0066CC', '99CCFF', '66B2FF', '66FFFF', '006633', '00CC66', '66FF66', '00FF00', '009900', 'FFFF99', 'FFFF00', 'CCCC00', 'FFB266', 'CC6600', '994C00', 'FF9999', 'FF0000', 'CC0000', '990000'}
myColorMap = zeros(length(hexMap), 3); % Preallocate
for k = 1 : length(hexMap)
	thisCell = hexMap{k}
	r = hex2dec(thisCell(1:2));
	g = hex2dec(thisCell(3:4));
	b = hex2dec(thisCell(5:6));
	myColorMap(k, :) = [r, g, b];
end
myColorMap = myColorMap / 255; % Normalize to range 0-1

colormap(myColorMap);
colorbar;
xticks([1:10])
yticks([1:19])

box off
%axis off
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 9 8]);
print(sprintf('%s/%s/fig/%s/%s/_map_val_zoom', directory_name, Fig_num,experiment_name,experiment_type), '-dsvg')

