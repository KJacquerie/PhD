clear all 
close all 
clc

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

%%

experiment_mat = ["Graupner2012"]%, "Graupner2012_RESET", "Graupner2012_noBURST"];  
directory_name = "/Users/kathleen/Documents/PhD/2023-Project/Fig2_MNIST"; 




%%

N_cycles = 16; 
NB_dgt = 10; 
tauG = 200; 
NB_samples=1; 
NB_states=N_cycles*2;
Duration_cycle = 30000; 

nPost = 10; 
nPre = 484; %144 
ncellsC = nPost + nPre;


dt=0.01;

NB_grid = 22;%12;
NB_pixels = NB_grid; 

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


mean_vec = zeros(NB_pixels^2,1); 

load('mean_MNIST.dat')
NB_pixels = sqrt(size(mean_MNIST,1)); 
for idx_dgt=1:1:NB_dgt
    X = reshape(mean_MNIST(:,idx_dgt), NB_pixels, NB_pixels); 
    mean_cell{idx_dgt} = (X-min(min(X)))/(max(max(X))-min(min(X))); 
    mean_vec((idx_dgt-1)*NB_pixels^2+1:idx_dgt*NB_pixels^2) = reshape(mean_cell{idx_dgt}, [],1)
end


%%

w_shaped  = {}; 
g_shaped  = {}; 
wg_shaped = {};
TONIC = {}; 
BURST = {}; 
clims = [0 1]; 
count=1; 

TYPE_NORM = 'REL'; % 'ABS';

wg_cell={}; 
chosen_STATE=[NB_states-1, NB_states]; 
for idx_expm = 1:1:length(experiment_mat)
    experiment_name = experiment_mat{idx_expm}; 
    w = load(sprintf('%s/data/%s/w_state.dat', directory_name,  experiment_name));
    g = load(sprintf('%s/data/%s/g_state.dat', directory_name,  experiment_name));
    dgt_presented = load(sprintf('%s/data/digit_presented.dat', directory_name));
    dgt_presented = dgt_presented-1; 
    
    wg_cell{idx_expm} = w.*g; 
    switch TYPE_NORM
        case 'REL'
            wgMAXX =max(max(w.*g));%0.0011; %;%1e-3;%4.3e-3;%max(max(w.*g))-0.001; 
            wgMIN = min(min(w.*g));%0;%9.6141e-06;%min(min(w.*g)); 
        case 'ABS'
            wgMAXX =0.0011; %;%1e-3;%4.3e-3;%max(max(w.*g))-0.001; 
            wgMIN = 0;%9.6141e-06;%min(min(w.*g)); 
    end
    
    
    count_fig=1; 
    for idx_chosen = 1:1:length(chosen_STATE)
        idx = chosen_STATE(idx_chosen); 
        figure(count_fig)
        count=1;
        for idx_dgt = 1:1:NB_dgt
            subplot(length(experiment_mat)+1, NB_dgt, count+(idx_expm-1)*10)
            w_shaped = reshape(w((idx_dgt-1)*nPre+1:idx_dgt*nPre,idx ),NB_grid,NB_grid);
            g_shaped = reshape(g((idx_dgt-1)*nPre+1:idx_dgt*nPre,idx ),NB_grid,NB_grid); 
            wg_shaped = (w_shaped.*g_shaped - wgMIN)/(wgMAXX - wgMIN);
            max(max(wg_shaped))
            imagesc(wg_shaped',clims)  
            colormap(gray(101))
            axis off
            box off 
            count=count+1;
            
        end
         
        count_fig=count_fig+1;
    end
    


end

for idx_dgt=1:1:NB_dgt
    figure(1)
    subplot(length(experiment_mat)+1, NB_dgt, count+(idx_expm-1)*10)
    imagesc(mean_cell{idx_dgt}',clims)  
    colormap(gray(101))
    axis off
    box off 
    
    figure(2)
    subplot(length(experiment_mat)+1, NB_dgt, count+(idx_expm-1)*10)
    imagesc(mean_cell{idx_dgt}',clims)  
    colormap(gray(101))
    axis off
    box off 
    
    count=count+1; 
end

set(figure(1),'PaperPositionMode','auto');
set(figure(1), 'PaperUnits', 'centimeters');
set(figure(1), 'PaperPosition', [0 0 50 5]);
%print(sprintf('%s/fig/LastTonic_RF',directory_name),'-depsc', '-painters')


set(figure(2),'PaperPositionMode','auto');
set(figure(2), 'PaperUnits', 'centimeters');
set(figure(2), 'PaperPosition', [0 0 50 5]);
%print(sprintf('%s/fig/LastBurst_RF',directory_name),'-depsc', '-painters')


%%
Nshown = 3*ones(1,NB_dgt); 

%for idx_dgt=1:1:NB_dgt
%    Nshown(idx_dgt)=length(find(seq==idx_dgt-1));
%end
%%



Err = {}; 
TEMP = zeros( NB_dgt, NB_states); 
for idx_expm=1:1:length(experiment_mat)
    for idx_state=1:1:NB_states
        for idx_dgt=1:1:NB_dgt
            
            %TEMP(idx_expm, idx_dgt) = mean((mean_cell{idx_dgt}-TONIC{idx_expm, idx_dgt}).^2, 'all');
              TEMP(idx_expm, idx_dgt) = immse(mean_cell{idx_dgt},TONIC{idx_expm, idx_dgt})
        end
    end
    Err{idx_state} = TEMP; 
end

