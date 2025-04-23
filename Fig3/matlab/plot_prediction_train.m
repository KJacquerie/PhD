clear all 
close all 
clc

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

%%

experiment_mat = [ "Graupner2012","Graupner2012_learn","Graupner2012_Reset", "Graupner2012_noBURST"];%, "Graupner2012_scale"]; 
%experiment_name = "Graupner2012"; 
directory_name = "/Users/kathleen/Documents/PhD/2023-Project/Fig3"; 
type_dataset = 'fair_trainingset'; 
NB_samples= 20;
NB_digits=4; 

%%


nPost = 10; 

color_blue = [106 153 208]./255; 
color_gray = [ 0.5 0.5 0.5]; 
color_pink = [250/255 244/255 247/255];
color_reset = [ 157 12 58]./255; 
color_corr = [63/255 92/255 206/255]; 
color_uncorr = [130/255 187/255 255/255]; 
color_green = [112 172 71]./255; 
color_mat = [200 52 93; 
            255 193 0;
            112 173 71;
            237 125 49
            91 155 213; 
            ]./255; 
        
        

%%
chosen_STATE=[1:100]; %%;%;%[1:100]; %;%[381:402]%[402-11:402]%[1:6]% 398 399 400 401 402]; 

%idx_presented = [1:20, 51:70, 101:120, 151:170, 201:220, 251:270, 301:320, 351:370, 401:420, 451:470]; %[1:1:500];


dgt_presented = zeros(1,NB_samples); 
for idx_loop=1:1:NB_digits-1
    dgt_presented = [dgt_presented idx_loop*ones(1,NB_samples)]; 
end


Prediction = zeros(length(experiment_mat), length(chosen_STATE));
for idx_expm = 1:1:length(experiment_mat)
    experiment_name = experiment_mat{idx_expm}; 
  
    count=1; 
    for idx_chosen = 1:1:length(chosen_STATE)
        %idx = chosen_STATE(idx_chosen); 
        %Spkt_mat = load(sprintf('%s/data/%s/%s/spkt_mat_%d.dat', directory_name, type_dataset, experiment_name, chosen_STATE(idx_chosen)));
        Spkt_mat = load(sprintf('%s/data/%s/spkt_mat_%d.dat', directory_name, experiment_name, chosen_STATE(idx_chosen)));
        Spkt{idx_expm}{idx_chosen} = Spkt_mat; 
        %figure(idx_expm)
        %subplot(1,length(chosen_STATE), idx_chosen)
        %imagesc(Spkt_mat, [0 90])
        %axis off
        for idx_dgt=1:1:size(Spkt_mat,2)
            
            [~,idx_MAX] = max(Spkt_mat(:,idx_dgt));
            Spkt_max{idx_expm}(idx_dgt,idx_chosen) = idx_MAX(1)-1;
            
            if((idx_MAX(1)-1)==dgt_presented(idx_dgt))
                Prediction(idx_expm, idx_chosen)= Prediction(idx_expm,idx_chosen)+1; 
            end
        end
    end
    %subplot(1,length(chosen_STATE)+1, idx_chosen+1)
    %axis off
        %olorbar
%         box off
%         grid on
%         ax=gca;
%         ax.XAxis.FontSize = 11;
%         ax.YAxis.FontSize = 11;
% 
%         set(gcf,'PaperPositionMode','auto');
%         set(gcf, 'PaperUnits', 'centimeters');
%         set(gcf, 'PaperPosition', [0 0 35 5]);
%         print(sprintf('%s/fig/%s_state%d', directory_name, experiment_name, chosen_STATE(idx_chosen)), '-dsvg')
%         
%    
end



%%
color_rule = ["347889", "78cabf","f68827", "da4300", "d8e150"]; %cyan-blue-orange-red-green
%color_rule = ["347889",  "9FB790", "f68827", "BF3520"]; % jaune "FEB215" % 347889 blue % green 679FB1

figure
hold on
for idx=1:1:size(Prediction,1)
    plot(Prediction(idx,:)*(100/size(Spkt_mat,2)), '-o','linewidth', 1.5, 'MarkerSize', 2,'color', hex2rgb(color_rule(idx)),'MarkerFaceColor',hex2rgb(color_rule(idx)), 'MarkerEdgeColor',hex2rgb(color_rule(idx)) )

end
ylim([0 100])
xlim([-0.2 size(Prediction,2)])
xticks([1:20:size(Prediction,2)])
xticklabels({''})
yticks([0 25 50 75 100])
yticklabels({''})
box off
grid on

ax=gca;
ax.XAxis.FontSize = 11;
ax.YAxis.FontSize = 11;
%ax.TickDir= 'out';



if(length(chosen_STATE)<40)
    set(gcf, 'PaperPosition', [0 0 3.4 5]); %17 - %3.7
    print(sprintf('%s/fig/%s/prediction_comp_end', directory_name, type_dataset), '-dsvg')

else
    set(gcf, 'PaperPosition', [0 0 16 5]); %17 - %3.7
    print(sprintf('%s/fig/%s/prediction_comp_start', directory_name, type_dataset), '-dsvg')

end



%%
