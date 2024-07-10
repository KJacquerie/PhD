clear all
close all
clc

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

directory_name = "/Users/kathleen/Documents/PhD/2023-Project"; 
experiment_name = "Graupner2012_VAR"; 
experiment_type = "low_dep"; 
bound_type = "SB";
fMax=50; 
new_curves=0; 
Fig_num = "Fig5"; 

tau_mat = 100:100:1000;

SNR_count=zeros(length(100:100:1000), 7)

%tau_val= 100;

for idx_tau = 1:1:size(tau_mat,2)
    tau_val = tau_mat(idx_tau)
for idx_lINIT=1:1:7
    idx_lINIT
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


    %% w
    %points = [1 19 32 47 80] ;
    figure
    subplot(4,1,1)
    hold on
    for idx_t=1:1:10
        plot([idx_t idx_t], [-1 1], 'linewidth', 0.25, 'color', [0.3 0.3 0.3]) 
    end
    plot(w_STATE','o', 'color', [0.75 0.75 0.75], 'linewidth', 0.25)
    for idx_p=1:1:5
        plot( w_STATE(idx_p,:)','o','color', 'k','linewidth', 1)
    end


    %xlim([T1 T2])
    mm = min(min(w_STATE)); 
    MM = max(max(w_STATE)); 
    ylim([mm MM]); 
    xlim([0 10])
    %xticks([0 15000 30000 45000 60000 75000 90000 105000 120000 135000 150000])
    %xticklabels({'','','','','','','','','','','',''})
    yticks([mm MM])
    %yticklabels({'', ''})


    %set(gcf,'PaperPositionMode','auto');
    %set(gcf, 'PaperUnits', 'centimeters');
    %set(gcf, 'PaperPosition', [0 0 7 4]);
    %print(sprintf('%s/%s/fig/%s/%s/w_%d_%d',directory_name, Fig_num, experiment_name, experiment_type, tau_val, idx_lINIT), '-dsvg', '-painters')


    %mm_MM = [mm MM]; 
    %csvwrite(sprintf('%s/%s/fig/tau=%d/mm_MM',directory_name, Fig_num, tau_val), mm_MM)

    %figure
    subplot(4,1,2)
    hold on
    for idx_t=1:1:10
        plot([idx_t idx_t], [-1 1], 'linewidth', 0.25, 'color', [0.3 0.3 0.3]) 
    end
    plot(l_STATE','o', 'color', [0.75 0.75 0.75], 'linewidth', 0.25)
    for idx_p=1:1:5
        plot( l_STATE(idx_p,:)','o','color', 'k','linewidth', 1)
    end


    %xlim([T1 T2])
    mm = min(min(l_STATE)); 
    MM = max(max(l_STATE)); 
    ylim([mm MM]); 
    xlim([0 10])
    %xticks([0 15000 30000 45000 60000 75000 90000 105000 120000 135000 150000])
    %xticklabels({'','','','','','','','','','','',''})
    yticks([mm MM])
    %yticklabels({'', ''})


    %set(gcf,'PaperPositionMode','auto');
    %set(gcf, 'PaperUnits', 'centimeters');
    %set(gcf, 'PaperPosition', [0 0 7 4]);
    %print(sprintf('%s/%s/fig/%s/%s/l_%d_%d',directory_name, Fig_num, experiment_name, experiment_type, tau_val, idx_lINIT), '-dsvg', '-painters')


    %figure
    subplot(4,1,3)
    hold on
    for idx_t=1:1:10
        plot([idx_t idx_t], [-1 1], 'linewidth', 0.25, 'color', [0.3 0.3 0.3]) 
    end
    plot(wl_STATE','o', 'color', [0.75 0.75 0.75], 'linewidth', 0.25)
    for idx_p=1:1:5
        plot( wl_STATE(idx_p,:)','o','color', 'k','linewidth', 1)
    end


    %xlim([T1 T2])
    mm = min(min(wl_STATE)); 
    MM = max(max(wl_STATE)); 
    ylim([mm MM]); 
    xlim([0 10])
    %xticks([0 15000 30000 45000 60000 75000 90000 105000 120000 135000 150000])
    %xticklabels({'','','','','','','','','','','',''})
    yticks([mm MM])
    %yticklabels({'', ''})


    %set(gcf,'PaperPositionMode','auto');
    %set(gcf, 'PaperUnits', 'centimeters');
    %set(gcf, 'PaperPosition', [0 0 7 4*5]);
    %print(sprintf('%s/%s/fig/%s/%s/_%d_%d',directory_name, Fig_num, experiment_name, experiment_type, tau_val, idx_lINIT), '-dsvg', '-painters')


    

    %%

    count=1; 
    for idx_state=1:2:Ncycles*2-1

        wl_tonic(:,count) = wl_STATE(:,idx_state);
        wl_burst(:,count) = wl_STATE(:,idx_state+1);

        SNR_tonic(count)  = max(wl_tonic(:,count))/mean(wl_tonic(:,count)); 
        SNR_burst(count) = max(wl_burst(:,count))/mean(wl_burst(:,count)); 
        
        if(SNR_tonic(count)>SNR_burst)
            SNR_count(idx_tau, idx_lINIT) = SNR_count(idx_tau, idx_lINIT) -1;
            SNR_param{idx_tau,idx_lINIT}(count)=-1;
        else
            if(SNR_tonic(count)<SNR_burst(count))
                SNR_count(idx_tau, idx_lINIT) = SNR_count(idx_tau, idx_lINIT) +1;
                SNR_param{idx_tau,idx_lINIT}(count)=+1;
            end
        end
        %SNR_tonic(count)  = mean(wl_tonic(1:5,count))/mean(wl_tonic(6:100,count)); 
        %SNR_burst(count) = mean(wl_burst(1:5,count))/mean(wl_burst(6:100,count)); 
        count=count+1; 
    end


    SNR = [SNR_tonic; SNR_burst];

    %figure
    %bar(SNR')
    subplot(4,1,4)
    hold on
    bar([1:2:2*size(SNR_tonic,2)], SNR_tonic(1:end),'FaceColor',color_tonic,'EdgeColor',color_tonic,'LineWidth',0.5, 'BarWidth', 0.25)
    bar([2:2:2*size(SNR_tonic,2)], SNR_burst(1:end),'FaceColor',color_burst,'EdgeColor',color_burst,'LineWidth',0.5, 'BarWidth', 0.25)
    xlim([0 11])
    box off
    %axis off
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 7 5*5]);
    print(sprintf('%s/%s/fig/%s/%s/SNR_%d_%d', directory_name, Fig_num,experiment_name,experiment_type,  tau_val, idx_lINIT), '-dsvg')


    close all
    
    %upselection 2
    %downselection -2
    % no change 0
    % from down to up 1
    % from up to down -1
    if(~isempty(sort([strfind(SNR_param{idx_tau,idx_lINIT}>=0, [0 1]) strfind(SNR_param{idx_tau,idx_lINIT}>=0, [1 0])])))
        if(SNR_param{idx_tau,idx_lINIT}(1)>0)
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
