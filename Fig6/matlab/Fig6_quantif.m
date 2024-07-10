clc
%clear all 
close all 
%%

ptx = 7; 
pty=5; 
N_cycles =8;
N_samples=1;
NB_dgt=2; 
NB_states=N_cycles*2; 

directory_name = "/Users/kathleen/Documents/PhD/2023-Project/Fig6"; 
model_mat = ["Graupner2012", "Graupner2016"]; 


experiment_mat = ["nOL", "OL"];



nPost = 1; 
nPre = 16; 
ncellsC = nPost + nPre;

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
dt=0.01;

NB_grid = 4; 

%%

clims = [0 1]; 
for idx_model=1:1:length(model_mat)
    for idx_expm=1:1:length(experiment_mat)
     model = model_mat(idx_model)
     experiment_name = experiment_mat(idx_expm)
    w = load(sprintf('%s/data/%s/%s/w_state%d.dat', directory_name, model, experiment_name, N_cycles));
    g = load(sprintf('%s/data/%s/%s/g_state%d.dat', directory_name, model, experiment_name, N_cycles));

    wg = w.*g;
    [wgMAXX, maxwg_index] = max(wg(:));
    [wgMIN, minwg_index] = min(wg(:));
    mm_MM = [min(w(:)) max(w(:))]; 
    csvwrite(sprintf('%s/fig/%s/%s/mm_MM',directory_name, model, experiment_name), mm_MM)
    
    l_mm_MM = [min(g(:)) max(g(:))]; 
    csvwrite(sprintf('%s/fig/%s/%s/l_mm_MM',directory_name, model, experiment_name), l_mm_MM)

    wl_mm_MM = [min(wg(:)) max(wg(:))]; 
    csvwrite(sprintf('%s/fig/%s/%s/wl_mm_MM',directory_name, model, experiment_name), wl_mm_MM)

    %%
    mean_cell={};
    switch experiment_name
        case "OL"
            mean_cell{1} = zeros(4,4); 
            mean_cell{2} = zeros(4,4);

            mean_cell{1}(:,1:2) = 1;
            mean_cell{2}(1:2,:) = 1;
        case "nOL"
            mean_cell{1} = zeros(4,4); 
            mean_cell{2} = zeros(4,4);

            mean_cell{1}(:,1:2) = 1;
            mean_cell{2}(:,3:4) = 1;
    end


%%


    for idx = 1:1:NB_states
        for idx_dgt = 1:1:NB_dgt
            w_shaped = reshape(w((idx_dgt-1)*nPre+1:idx_dgt*nPre,idx ),NB_grid,NB_grid);
            g_shaped = reshape(g((idx_dgt-1)*nPre+1:idx_dgt*nPre,idx ),NB_grid,NB_grid); 
            wg_shaped = (w_shaped.*g_shaped - wgMIN)/(wgMAXX - wgMIN);
            %Temp = immse(mean_cell{idx_dgt}, w_shaped.*g_shaped); 
            
            %Temp = immse(mean_cell{idx_dgt},wg_shaped); 
            TempC = corr2(mean_cell{idx_dgt},wg_shaped); 
            
            %Err_Model(idx, idx_dgt) = Temp; 
            Corr_Model(idx, idx_dgt) = TempC; 
            
        end
         
    end

switch experiment_name
    case "nOL"
        idxx = 1;
    case "OL"
        idxx = 2;
end

switch model
    case "Graupner2012"
        idxy=1; 
    case "Graupner2016"
        idxy=2; 
end

%Err_cell{idxx, idxy}= Err_Model;
Corr_cell{idxx, idxy}= Corr_Model;
%csvwrite(sprintf('%s/fig/Err_Model',directory_name), Err_Model)
%csvwrite(sprintf('%s/fig/Corr_Model',directory_name), Corr_Model)
  
    end
end
%%


% for idx_dgt=1:1:NB_dgt
%     figure(idx_dgt)
%     hold on
%     plot([1:1:NB_states], Err_cell{1,1}(:,idx_dgt), '-o', 'linewidth',3)%, 'color', [0 0 1])
%     plot([1:1:NB_states], Err_cell{1,2}(:,idx_dgt), '-o', 'linewidth',3)%, 'color', 'r')%[1 0.5 0])
%     plot([1:1:NB_states], Err_cell{2,1}(:,idx_dgt), ':o', 'linewidth',3)%, 'color', 'g')%[0.7 0 1])
%     plot([1:1:NB_states], Err_cell{2,2}(:,idx_dgt), ':o', 'linewidth',3)%, 'color', [0.5 0.5 0.5])
%     
%     xlim([0 17])
%     set(figure(idx_dgt),'PaperPositionMode','auto');
%     set(figure(idx_dgt), 'PaperUnits', 'centimeters');
%     set(figure(idx_dgt), 'PaperPosition', [0 0 30 15]);
%     xlabel('state')
%     ylabel('MSE')
%     print(sprintf('%s/fig/ERR_%d',directory_name, idx_dgt),'-dsvg', '-painters')
% end

%%

% 
% for idx_dgt=1:1:NB_dgt
%     figure(idx_dgt)
%     hold on
%     plot([1:1:NB_states], Corr_cell{1,1}(:,idx_dgt), '-o', 'linewidth',3)%, 'color', [0 0 1])
%     plot([1:1:NB_states], Corr_cell{1,2}(:,idx_dgt), '-o', 'linewidth',3)%, 'color', 'r')%[1 0.5 0])
%     plot([1:1:NB_states], Corr_cell{2,1}(:,idx_dgt), ':o', 'linewidth',3)%, 'color', 'g')%[0.7 0 1])
%     plot([1:1:NB_states], Corr_cell{2,2}(:,idx_dgt), ':o', 'linewidth',3)%, 'color', [0.5 0.5 0.5])
%     
%     xlim([0 17])
%     set(figure(idx_dgt),'PaperPositionMode','auto');
%     set(figure(idx_dgt), 'PaperUnits', 'centimeters');
%     set(figure(idx_dgt), 'PaperPosition', [0 0 30 15]);
%     xlabel('state')
%     ylabel('Corr')
%     print(sprintf('%s/fig/Corr_%d',directory_name, idx_dgt),'-dsvg', '-painters')
% end


%%
color = [ "E7A7A1","3C185A"]'; 

count=1; 
for idx_expm=1:1:length(experiment_mat)
    for idx_model=1:1:length(model_mat)
        figure(count)
        hold on
        for idx_dgt=1:1:2
            plot([1:1:NB_states], Corr_cell{idx_expm,idx_model}(:,idx_dgt), 'color', hex2rgb(color(idx_dgt,:)), 'linewidth', 1.5); 
        end
        
        ylim([-0.6 1.1])
        xlim([0 17])
        set(figure(count),'PaperPositionMode','auto');
        set(figure(count), 'PaperUnits', 'centimeters');
        set(figure(count), 'PaperPosition', [0 0 8 4]);
        %xticks([1:16])
        %xticklabels({'','','','','','','','','','','','','','','',''})
        xticks([1:8:16])
        xticklabels({'','','','','','','','','','','','','','','',''})
        yticks([-0.5:0.5:1])
        yticklabels({'','', '', ''})
        print(sprintf('%s/fig/Corr_%s_%s',directory_name, model_mat(idx_model),experiment_mat(idx_expm)),'-dsvg', '-painters')
        count=count+1;
    end
end

%%

%gray_sq = rand(NB_grid, NB_grid);
% Err_rand = zeros(NB_dgt,1);
% Corr_rand = zeros(NB_dgt,1);
% 
% close all
% for idx_dgt=1:1:NB_dgt
%     Err_rand(idx_dgt) = immse(mean_cell{idx_dgt}, gray_sq); 
%     Corr_rand(idx_dgt) = corr2(mean_cell{idx_dgt}, gray_sq); 
%     figure(idx_dgt)
%     bar([1,2], [Err_rand(idx_dgt), Corr_rand(idx_dgt)])
%     set(figure(idx_dgt),'PaperPositionMode','auto');
%     set(figure(idx_dgt), 'PaperUnits', 'centimeters');
%     set(figure(idx_dgt), 'PaperPosition', [0 0 5 5]);
%     ylim([-0.15 0.45])
%     %xlabel('state')
%     %ylabel('MSE')
%     print(sprintf('%s/fig/RAND_%d',directory_name, idx_dgt),'-dsvg', '-painters')
% end

%%

figure
imagesc(mean_cell{1,2})
colormap(gray)
colorbar