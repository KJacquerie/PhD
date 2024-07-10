clear all 
close all 
clc

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

%%

experiment_mat = ["Graupner2012", "Graupner2012_learn", "Graupner2012_RESET", "Graupner2012_noBURST"]; 
%experiment_name = "Graupner2012"; 
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

NB_grid = 22;%12;

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

%%

w_shaped  = {}; 
g_shaped  = {}; 
wg_shaped = {};
TONIC = {}; 
BURST = {}; 
clims = [0 1]; 
count=1; 

TYPE_NORM = 'REL'; % 'ABS' / 'REL';

chosen_STATE=[NB_states-3, NB_states-2, NB_states-1, NB_states]; 
for idx_expm = 1:1:length(experiment_mat)
    experiment_name = experiment_mat{idx_expm}; 
    w = load(sprintf('%s/data/%s/w_state.dat', directory_name,  experiment_name));
    g = load(sprintf('%s/data/%s/g_state.dat', directory_name,  experiment_name));
    dgt_presented = load(sprintf('%s/data/digit_presented.dat', directory_name));
    dgt_presented = dgt_presented-1; 
    
    switch TYPE_NORM
        case 'REL'
            wgMAXX =max(max(w.*g));%0.0011; %;%1e-3;%4.3e-3;%max(max(w.*g))-0.001; 
            wgMIN = min(min(w.*g));%0;%9.6141e-06;%min(min(w.*g)); 
        case 'ABS'
            wgMAXX =0.001; %;%1e-3;%4.3e-3;%max(max(w.*g))-0.001; 
            wgMIN = 0;%9.6141e-06;%min(min(w.*g)); 
    end
    
    
    count_fig=1; 
    for idx_chosen = 1:1:length(chosen_STATE)
        idx = chosen_STATE(idx_chosen); 
        figure(count_fig)
        count=1;
        for idx_dgt = 1:1:NB_dgt
            %wgMAXX =max(max(w((idx_dgt-1)*nPre+1:idx_dgt*nPre,idx ).*g((idx_dgt-1)*nPre+1:idx_dgt*nPre,idx )));
            %wgMIN = min(min(w((idx_dgt-1)*nPre+1:idx_dgt*nPre,idx ).*g((idx_dgt-1)*nPre+1:idx_dgt*nPre,idx )))
            subplot(length(experiment_mat)+1, NB_dgt, count+(idx_expm-1)*10)
            w_shaped = reshape(w((idx_dgt-1)*nPre+1:idx_dgt*nPre,idx ),NB_grid,NB_grid);
            g_shaped = reshape(g((idx_dgt-1)*nPre+1:idx_dgt*nPre,idx ),NB_grid,NB_grid); 
            wg_shaped = (w_shaped.*g_shaped - wgMIN)/(wgMAXX - wgMIN);
            max(max(wg_shaped))
            %if(mod(idx, 2)==0)
            %    BURST{idx_expm, idx_dgt} = wg_shaped; 
            %else
            %    TONIC{idx_expm, idx_dgt} = wg_shaped; 
            %end
            imagesc(wg_shaped',clims)  
            colormap(gray(101))
            axis off
            box off 
            count=count+1;
            
        end
         
        count_fig=count_fig+1;
    end
   
end

for idx_fig=1:1:length(chosen_STATE)
    set(figure(idx_fig),'PaperPositionMode','auto');
    set(figure(idx_fig), 'PaperUnits', 'centimeters');
    set(figure(idx_fig), 'PaperPosition', [0 0 50 5]);
    print(sprintf('%s/fig/LastTonic_RF',directory_name),'-depsc', '-painters')
end


for idx_fig=1:1:length(chosen_STATE)
    count=NB_dgt+1; 
    for idx_dgt=1:1:NB_dgt
        figure(idx_fig)
        subplot(length(experiment_mat)+1, NB_dgt,count+(idx_expm-1)*10)
        imagesc(mean_cell{idx_dgt}',clims)  
        colormap(gray(101))
        axis off
        box off 

        count=count+1; 
    end
end

set(figure(1),'PaperPositionMode','auto');
set(figure(1), 'PaperUnits', 'centimeters');
set(figure(1), 'PaperPosition', [0 0 50 30]);
print(sprintf('%s/fig/LastTonic_RF',directory_name),'-dsvg', '-painters')


set(figure(2),'PaperPositionMode','auto');
set(figure(2), 'PaperUnits', 'centimeters');
set(figure(2), 'PaperPosition', [0 0 50 30]);
print(sprintf('%s/fig/LastBurst_RF',directory_name),'-dsvg', '-painters')



%%
% Nshown = 3*ones(1,NB_dgt);%zeros(1,NB_dgt); 
% seq = dgt_presented(1:2:end); 
% for idx_dgt=1:1:NB_dgt
%     Nshown(idx_dgt)=length(find(seq==idx_dgt-1));
% end
%%
% Err_TONIC = zeros(length(experiment_mat), NB_dgt); 
% for idx_expm=1:1:size(TONIC,1)
%     for idx_dgt=1:1:NB_dgt
%         Err_TONIC(idx_expm, idx_dgt) = mean((mean_cell{idx_dgt}-TONIC{idx_expm, idx_dgt}).^2, 'all');
%     end
% end
% 
% Err_BURST = zeros(length(experiment_mat), NB_dgt); 
% for idx_expm=1:1:size(TONIC,1)
%     for idx_dgt=1:1:NB_dgt
%         Err_BURST(idx_expm, idx_dgt) = mean((mean_cell{idx_dgt}-BURST{idx_expm, idx_dgt}).^2, 'all');
%     end
% end
% 
% figure
% hold
% plot([0:1:9], Err_TONIC','-o')
% plot([0:1:9], Err_BURST','-x')
% 

close all 
%% 
Err_Model1=zeros(NB_states, NB_dgt);
Err_Model2=zeros(NB_states, NB_dgt);
Err_Model3=zeros(NB_states, NB_dgt);
Err_Model4=zeros(NB_states, NB_dgt);


color = ["F2CEB0", "E7A7A1","DE8293", "D55E86", "AE4B88", "88398D", "602892", "4E1F76", "3C185A", "29113E"]'; 
color_blue = ["487E91",  "679FB1","95D1E2", "9FB790", "BFDCAE", "F1D350","FEB215", "EF857D", "E94125","BF3520"];  

color_model = ["679FB1",  "9FB790", "FEB215", "BF3520"];
for idx_expm=1:1:length(experiment_mat)
    experiment_name = experiment_mat{idx_expm}; 
    w = load(sprintf('%s/data/%s/w_state.dat', directory_name,  experiment_name));
    g = load(sprintf('%s/data/%s/g_state.dat', directory_name,  experiment_name));
    wgMAXX = max(max(w.*g));%0.0011; %;%1e-3;%4.3e-3;%max(max(w.*g))-0.001; 
    wgMIN = min(min(w.*g));
    for idx = 1:1:NB_states
        for idx_dgt = 1:1:NB_dgt
            w_shaped = reshape(w((idx_dgt-1)*nPre+1:idx_dgt*nPre,idx ),NB_grid,NB_grid);
            g_shaped = reshape(g((idx_dgt-1)*nPre+1:idx_dgt*nPre,idx ),NB_grid,NB_grid); 
            wg_shaped = (w_shaped.*g_shaped - wgMIN)/(wgMAXX - wgMIN);
            %Temp = immse(mean_cell{idx_dgt}, w_shaped.*g_shaped); 
            %Temp = immse(mean_cell{idx_dgt},wg_shaped); 
            
            
            if(idx_expm==2)
                Temp = corr2(mean_learn{idx_dgt},wg_shaped); 
            else
                Temp = corr2(mean_cell{idx_dgt},wg_shaped); 
            end

            switch idx_expm
                case 1
                    Err_Model1(idx, idx_dgt) = Temp; 
                case 2
                    Err_Model2(idx, idx_dgt) = Temp; 
                case 3
                    Err_Model3(idx, idx_dgt) = Temp; 
                case 4
                    Err_Model4(idx, idx_dgt) = Temp; 
            end
            
        end
         
    end

end
   
 
for idx_dgt=1:1:NB_dgt
    figure(idx_dgt)
    hold on

  %  if(dgt_presented(idx_dgt)==idx_dgt && mod(idx_dgt,2)==1)
%         plot([1 NB_states], [0 0], '-', 'color', [0.7 0.7 0.7])
%     end
%    plot([1 NB_states], [0 0], '-', 'color', [0.7 0.7 0.7])
    plot([1:1:NB_states], Err_Model1(:,idx_dgt), '-o', 'color', hex2rgb(color_model(1)), 'MarkerFaceColor', hex2rgb(color_model(1)))
    plot([1:1:NB_states], Err_Model2(:,idx_dgt), '-o', 'color', hex2rgb(color_model(2)), 'MarkerFaceColor', hex2rgb(color_model(2)))
    plot([1:1:NB_states], Err_Model3(:,idx_dgt), '-o', 'color', hex2rgb(color_model(3)), 'MarkerFaceColor', hex2rgb(color_model(3)))
    plot([1:1:NB_states], Err_Model4(:,idx_dgt), '-o', 'color', hex2rgb(color_model(4)), 'MarkerFaceColor', hex2rgb(color_model(4)))
    
    x = find(dgt_presented==idx_dgt-1);
    if(isnan(x)==0)
        for idx_x=1:1:length(x)
            if(mod(idx_x,2)==1)
                plot([x(idx_x) x(idx_x)], [0 .5], '-', 'color', [0.5 0.5 0.5])
            end
        end
    end
    xlim([0 402])
    set(figure(idx_dgt),'PaperPositionMode','auto');
    set(figure(idx_dgt), 'PaperUnits', 'centimeters');
    set(figure(idx_dgt), 'PaperPosition', [0 0 30 15]);
    xlabel('state')
    ylabel('CORR')
    print(sprintf('%s/fig/_CORR_%d',directory_name, idx_dgt),'-dsvg', '-painters')
end

%%
for idx_model=1:1:4
    switch idx_model
        case 1
            Err_Mod = Err_Model1;
        case 2 
            Err_Mod = Err_Model2;
        case 3
            Err_Mod = Err_Model3;
        case 4
            Err_Mod = Err_Model4;
    end
    figure(100+idx_model)
    cmap = colormap(parula(10));
    hold on
    for idx_dgt=1:1:NB_dgt  
        plot(1:402, Err_Mod(:,idx_dgt), 'linewidth', 1.2, 'color', hex2rgb(color(idx_dgt))) %, 'color', cmap(idx_dgt,:),
    end
    yticks(-1:1:1)
    yticklabels({'', '', '', '', ''})
    xticks(1:20:402)
    xticklabels({'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ''})
    set(figure(100+idx_model),'PaperPositionMode','auto');
    set(figure(100+idx_model), 'PaperUnits', 'centimeters');
    set(figure(100+idx_model), 'PaperPosition', [0 0 12 5]);
    xlim([1 410])
    ylim([-1.1 1.1])
    %xlabel('state')
    %ylabel('CORR')
    print(sprintf('%s/fig/_CORR_Model_%d',directory_name, idx_model),'-dsvg', '-painters')
    
    xticks(1:1:22)
    xticklabels({'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ''})
    set(figure(100+idx_model),'PaperPositionMode','auto');
    set(figure(100+idx_model), 'PaperUnits', 'centimeters');
    set(figure(100+idx_model), 'PaperPosition', [0 0 4 5]);
    xlim([1 20])
    ylim([-1.1 1.1])
    %xlabel('state')
    %ylabel('CORR')
    print(sprintf('%s/fig/_CORR_Model_%d_start',directory_name, idx_model),'-dsvg', '-painters')
    
    xticks(381:1:402)
    xticklabels({'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ''})
    
    set(figure(100+idx_model),'PaperPositionMode','auto');
    set(figure(100+idx_model), 'PaperUnits', 'centimeters');
    set(figure(100+idx_model), 'PaperPosition', [0 0 4 5]);
    xlim([381 402])
    ylim([-1.1 1.1])
    %xlabel('state')
    %ylabel('CORR')
    print(sprintf('%s/fig/_CORR_Model_%d_end',directory_name, idx_model),'-dsvg', '-painters')
    
    
end
%%
gray_sq = rand(NB_pixels, NB_pixels);
Err_rand = zeros(NB_dgt,1);
Corr_rand = zeros(NB_dgt,1);

imagesc(reshape(gray_sq,NB_grid, NB_grid))

%%
 
% close all
% for idx_dgt=1:1:NB_dgt
%     Err_rand(idx_dgt) = immse(mean_cell{idx_dgt}, gray_sq); 
%     Corr_rand(idx_dgt) = corr2(mean_cell{idx_dgt}, gray_sq); 
%     figure(idx_dgt)
%     bar([1,2], [Err_rand(idx_dgt), Corr_rand(idx_dgt)])
%     set(figure(idx_dgt),'PaperPositionMode','auto');
%     set(figure(idx_dgt), 'PaperUnits', 'centimeters');
%     set(figure(idx_dgt), 'PaperPosition', [0 0 5 8]);
%     ylim([-0.05 0.35])
%     %xlabel('state')
%     %ylabel('MSE')
%     print(sprintf('%s/fig/_RAND_%d',directory_name, idx_dgt),'-dsvg', '-painters')
% end

%% correlation entre les nombres

Corr_dgt = zeros(NB_dgt, NB_dgt); 

figure
for idx_dgt=1:1:NB_dgt
    for idx_2=1:1:NB_dgt
        Corr_dgt(idx_dgt, idx_2)=corr2(mean_cell{idx_dgt}, mean_cell{idx_2}); 

    end
end

imagesc(Corr_dgt)
xticks=[1:1:10];
xticklabels({'0','1','2','3','4','5','6','7','8','9'})
yticks=[1:1:10];
yticklabels({'0','1','2','3','4','5','6','7','8','9'})
colorbar