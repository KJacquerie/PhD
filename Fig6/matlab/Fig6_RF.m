clear all 
close all 
clc

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');


%%

print_wg =1; 
pt= 11; 
ptx = 7; 
pty=5; 
N_cycles =8;
N_samples=1;
NB_dgt=2; 
NB_states=N_cycles*2; 

directory_name = "/Users/kathleen/Documents/PhD/2023-Project/Fig6"; 
model = "Graupner2016";
experiment_name =  "OL"; % nOL


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

w = load(sprintf('%s/data/%s/%s/w_state%d.dat', directory_name, model, experiment_name, N_cycles));
g = load(sprintf('%s/data/%s/%s/g_state%d.dat', directory_name, model, experiment_name, N_cycles));

wg = w.*g;
[max_wg, maxwg_index] = max(wg(:));
[min_wg, minwg_index] = min(wg(:));

[max_w, maxw_index] = max(w(:));
[min_w, minw_index] = min(w(:));
    


count=1;
xx=NB_states; 

figure(1) 
tiledlayout(1,xx, 'TileSpacing', 'non', 'Padding', 'none')

figure(2) 
tiledlayout(1,xx, 'TileSpacing', 'non', 'Padding', 'none')
wgMIN  = min(min(w.*g));
wgMAXX = max(max(w.*g));

mm_MM = [wgMIN wgMAXX]; 
csvwrite(sprintf('%s/fig/%s/%s/mm_MM',directory_name, model, experiment_name), mm_MM)
     
%%
    for idx = 1:1:xx%NB_states
            %disp(idx_state)
            
            w1 = reshape(w(1:nPre,idx),NB_grid,NB_grid); 
            g1 = reshape(g(1:nPre,idx),NB_grid,NB_grid); 
            
            w2 = reshape(w((nPre+1):2*nPre, idx),NB_grid,NB_grid);
            g2 = reshape(g((nPre+1):2*nPre, idx),NB_grid,NB_grid);
            w1_reset = (w1 -wgMIN)/(wgMAXX-wgMIN);
            wg1 = (w1.*g1 - wgMIN)/(wgMAXX - wgMIN);
            wg2 = (w2.*g2 - wgMIN)/(wgMAXX - wgMIN);
 
             
            figure(1)
            nexttile
            imagesc(wg1,clims)  
            colormap(gray(101))
            axis off
            box off

            figure(2)
            nexttile
            imagesc(wg2,clims)  
            colormap(gray(101))
            axis off
            box off
            
            
            ax=gca;
            ax.Position = 1.1*ax.Position; 
            count=count+1;
            
    end    
    set(figure(1),'PaperPositionMode','auto');
    set(figure(1), 'PaperUnits', 'centimeters');
    set(figure(1), 'PaperPosition', [0 0 50 5]);
    print(sprintf('%s/fig/%s/%s/RF1',directory_name, model, experiment_name),'-dsvg', '-painters')
           
    set(figure(2),'PaperPositionMode','auto');
    set(figure(2), 'PaperUnits', 'centimeters');
    set(figure(2), 'PaperPosition', [0 0 50 5]);
    print(sprintf('%s/fig/%s/%s/RF2',directory_name, model, experiment_name),'-dsvg', '-painters')
           


%% RF for each digit at each state separetly

close all 

clims = [0 1]; 
w_shaped  = {}; 
g_shaped  = {}; 
wg_shaped = {};


for idx_dgt = 1:1:NB_dgt
    
    for idx = 1:1:NB_states
        w_shaped{idx_dgt,idx} = reshape(w((idx_dgt-1)*nPre+1:idx_dgt*nPre,idx ),NB_grid,NB_grid);
        g_shaped{idx_dgt,idx} = reshape(g((idx_dgt-1)*nPre+1:idx_dgt*nPre,idx ),NB_grid,NB_grid); 
        wg_shaped{idx_dgt,idx} = (w_shaped{idx_dgt,idx}.*g_shaped{idx_dgt,idx} - wgMIN)/(wgMAXX - wgMIN);

        
        figure
        imagesc(wg_shaped{idx_dgt,idx},clims)  
        colormap(gray(101))
        axis off
        box off 
        set(gcf,'PaperPositionMode','auto');
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 1.5 1.5]);
        print(sprintf('%s/fig/%s/%s/RF%d_%d',directory_name, model, experiment_name, idx_dgt, idx),'-dsvg', '-painters')
        close all
    end
    

end

%%

chosen_state = [1:16];


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

        
        %nexttile
        imagesc(wg_shaped{idx_dgt,idx},clims)  
        colormap(gray(101))
        axis off
        box off 
        gcff = set(gcf); 
        gcff.WindowState = 'maximized';
        count=count+1;
    end
    
end
set(figure(100),'PaperPositionMode','auto');
set(figure(100), 'PaperUnits', 'centimeters');
set(figure(100), 'PaperPosition', [0 0 20 6]);
print(sprintf('%s/fig/%s/%s/STATE_RF',directory_name, model, experiment_name),'-dpdf', '-painters')
print(sprintf('%s/fig/%s/%s/STATE_RF',directory_name, model, experiment_name),'-dsvg', '-painters')



%% quantification pixel 

switch experiment_name
    case "nOL"  
        QUANTIF_pat = [mean(wg_shaped{1,end}(:,1:2),'all') mean(wg_shaped{1,end}(:,3:4),'all'); mean(wg_shaped{2,end}(:,1:2),'all') mean(wg_shaped{2,end}(:,3:4),'all')]
        %perc_pat1 = (QUANTIF_pat(1,1)- QUANTIF_pat(1,2))/QUANTIF_pat(1,2)
        %perc_pat2 = (QUANTIF_pat(2,2)- QUANTIF_pat(2,1))/QUANTIF_pat(2,1)
        
        perc_pat1 = 2*(QUANTIF_pat(1,1)- QUANTIF_pat(1,2))/(QUANTIF_pat(1,2)+QUANTIF_pat(1,1))
        perc_pat2 = 2*(QUANTIF_pat(2,1)- QUANTIF_pat(2,2))/(QUANTIF_pat(2,2)+QUANTIF_pat(2,1))
       
    case "OL"
        QUANTIF_pat = [mean(wg_shaped{1,end}(:,1:2),'all') mean(wg_shaped{1,end}(:,3:4),'all'); mean(wg_shaped{2,end}(1:2,:),'all') mean(wg_shaped{2,end}(3:4,:),'all')]
        %perc_pat1 = (QUANTIF_pat(1,1)- QUANTIF_pat(1,2))/QUANTIF_pat(1,2) 
        %perc_pat2 = (QUANTIF_pat(2,1)- QUANTIF_pat(2,2))/QUANTIF_pat(2,2) 
        
        perc_pat1 = 2*(QUANTIF_pat(1,1)- QUANTIF_pat(1,2))/(QUANTIF_pat(1,2)+QUANTIF_pat(1,1))
        perc_pat2 = 2*(QUANTIF_pat(2,1)- QUANTIF_pat(2,2))/(QUANTIF_pat(2,2)+QUANTIF_pat(2,1))
        
end

csvwrite(sprintf('%s/fig/%s/%s/perc_pat',directory_name, model, experiment_name),[perc_pat1 perc_pat2])
   
