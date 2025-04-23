clear all 
close all 
clc

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

%%


experiment_name = "Graupner2012_silent_1Hz"; 
directory_name = "/Users/kathleen/Documents/PhD/2023-Project/Fig2_MNIST_long"; 




%%

N_cycles = 201; 
NB_dgt = 10; 
tauG = 400; 
NB_samples=1; 
NB_states=N_cycles*2;
Duration_cycle = 30000; 

nPost = 10; 
nPre = 484; %144
ncellsC = nPost + nPre;


dt=0.01;

NB_grid = 22;%12

Duration_state = Duration_cycle/2; 
Duration_sleep = Duration_state; 
Duration_sample = Duration_state/NB_samples; 
Tdt_cycle = Duration_sleep + Duration_state; 
Tdt_state = Duration_state/dt;


print_wg =1; 
pt= 11; 
ptx = 7; 
pty=5; 


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

%%  SUBPLOT


clims = [0 1]; 

close all 

w = load(sprintf('%s/data/%s/w_state.dat', directory_name,  experiment_name));
g = load(sprintf('%s/data/%s/g_state.dat', directory_name,  experiment_name));
dgt_presented = load(sprintf('%s/data/digit_presented.dat', directory_name));


wgMAXX = max(max(w.*g)); %1e-3;%4.3e-3;%max(max(w.*g))-0.001; 
wgMIN  = min(min(w.*g)); 

mm_MM = [wgMIN wgMAXX]; 
csvwrite(sprintf('%s/fig/%s/mm_MM',directory_name, experiment_name), mm_MM)


xx=NB_states;
count=1; 
for idx=1:1:NB_dgt
    figure(idx) 
    tiledlayout(1,xx, 'TileSpacing', 'non', 'Padding', 'none')
end
close all

%% RF for each digt - subplot for each state

% w_shaped  = {}; 
% g_shaped  = {}; 
% wg_shaped = {};
%   
% for idx_dgt = 1:1:NB_dgt
%     figure(idx_dgt)
%     for idx = 1:1:NB_states
%         w_shaped{idx_dgt,idx} = reshape(w((idx_dgt-1)*nPre+1:idx_dgt*nPre,idx ),NB_grid,NB_grid);
%         g_shaped{idx_dgt,idx} = reshape(g((idx_dgt-1)*nPre+1:idx_dgt*nPre,idx ),NB_grid,NB_grid); 
%         wg_shaped{idx_dgt,idx} = (w_shaped{idx_dgt,idx}.*g_shaped{idx_dgt,idx} - wgMIN)/(wgMAXX - wgMIN);
% 
%         
%         nexttile
%         imagesc(wg_shaped{idx_dgt,idx}',clims)  
%         colormap(gray(101))
%         axis off
%         box off 
%     end
%     set(figure(idx_dgt),'PaperPositionMode','auto');
%     set(figure(idx_dgt), 'PaperUnits', 'centimeters');
%     set(figure(idx_dgt), 'PaperPosition', [0 0 50 5]);
%     print(sprintf('%s/fig/%s/RF%d',directory_name, experiment_name, idx_dgt),'-dsvg', '-painters')
% 
% end


%% RF for each digit at each state separetly

chosen_state = [1 2 201 202 219 220 381 382 399 400 401 402];
%chosen_state = [1:20:160]
  
for idx_dgt = 10:1:NB_dgt
    
    for i = 1:1:length(chosen_state)
        idx=chosen_state(i); 
        w_shaped{idx_dgt,idx} = reshape(w((idx_dgt-1)*nPre+1:idx_dgt*nPre,idx ),NB_grid,NB_grid);
        g_shaped{idx_dgt,idx} = reshape(g((idx_dgt-1)*nPre+1:idx_dgt*nPre,idx ),NB_grid,NB_grid); 
        wg_shaped{idx_dgt,idx} = (w_shaped{idx_dgt,idx}.*g_shaped{idx_dgt,idx} - wgMIN)/(wgMAXX - wgMIN);

        figure
        imagesc(wg_shaped{idx_dgt,idx}',clims)  
        colormap(gray(101))
        axis off
        box off 
        set(gcf,'PaperPositionMode','auto');
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 3 3]);
        print(sprintf('%s/fig/%s/RF%d_%d',directory_name, experiment_name, idx_dgt, idx),'-dsvg', '-painters')
        close all
    end
    

end

%%


load('mean_MNIST_20.dat')
NB_pixels = sqrt(size(mean_MNIST_20,1)); 
for idx_dgt=1:1:NB_dgt
    X = reshape(mean_MNIST_20(:,idx_dgt), NB_pixels, NB_pixels); 
    mean_cell{idx_dgt} = (X-min(min(X)))/(max(max(X))-min(min(X))); 
end


load('mean_MNIST_learn.dat')
NB_pixels = sqrt(size(mean_MNIST_learn,1)); 
for idx_dgt=1:1:NB_dgt
    X = reshape(mean_MNIST_learn(:,idx_dgt), NB_pixels, NB_pixels); 
    mean_learn{idx_dgt} = (X-min(min(X)))/(max(max(X))-min(min(X))); 
end

for idx_dgt=1:1:NB_dgt
    figure
    imagesc(mean_cell{idx_dgt}',clims)  
    colormap(gray(101))
    axis off
    box off 
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 3 3]);
    print(sprintf('%s/fig/RF_mean_%d',directory_name, idx_dgt),'-dsvg', '-painters')

    
    figure
    imagesc(mean_learn{idx_dgt}',clims)  
    colormap(gray(101))
    axis off
    box off 
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 3 3]);
    print(sprintf('%s/fig/RF_mean_learn_%d',directory_name, idx_dgt),'-dsvg', '-painters')
    close all
end

%% RF for each digt - subplot for each state

%chosen_state = [1 2 19 20 21 22 39 40 41 42 59 60 61 62];


chosen_state = [1 2 201 202 219 220 381 382 399 400 401 402];
chosen_state = [1:50:260]
w_shaped  = {}; 
g_shaped  = {}; 
wg_shaped = {};
  
figure(100)

count=1; 
for idx_dgt = 1:1:NB_dgt
    for idx_chosen = 1:1:length(chosen_state)
        subplot(NB_dgt, length(chosen_state),count)
        idx=chosen_state(idx_chosen);
        w_shaped{idx_dgt,idx} = reshape(w((idx_dgt-1)*nPre+1:idx_dgt*nPre,idx ),NB_grid,NB_grid);
        g_shaped{idx_dgt,idx} = reshape(g((idx_dgt-1)*nPre+1:idx_dgt*nPre,idx ),NB_grid,NB_grid); 
        wg_shaped{idx_dgt,idx} = (w_shaped{idx_dgt,idx}.*g_shaped{idx_dgt,idx} - wgMIN)/(wgMAXX - wgMIN);

        % Add a small constant to avoid log(0) issues
        %wg_shaped{idx_dgt, idx}(wg_shaped{idx_dgt, idx} == 0) = 0.0001;
        
        %nexttile
        imagesc(wg_shaped{idx_dgt,idx}',clims)  
        colormap(gray(101))
        axis off
        box off 
        gcff = set(gcf); 
        gcff.WindowState = 'maximized';
         % Set logarithmic color scale
        set(gca, 'ColorScale', 'log')
        
        count=count+1;
    end
    
end
set(figure(100),'PaperPositionMode','auto');
set(figure(100), 'PaperUnits', 'centimeters');
set(figure(100), 'PaperPosition', [0 0 50 55]);
print(sprintf('%s/fig/%s/STATE_RF',directory_name, experiment_name),'-dsvg', '-painters')







