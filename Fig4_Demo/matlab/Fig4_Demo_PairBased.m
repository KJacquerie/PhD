clear all
close all
clc

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

directory_name = "/Users/kathleen/Documents/PhD/2023-reset"; 
experiment_name = "PairBased"; 
bound_type = "SB";
region = "HPC"; 
Fig_num= "Fig4_Demo"; 
new_curves=1; 

Ncycles=1; 
Tcycle=100000; 
Tstate=50000; 
dt=0.01;
T = Ncycles*Tcycle;
Tdt = T/dt;
t = dt:dt:T;

tmp =  load(sprintf('%s/%s/data/idx_w.mat', directory_name, Fig_num));
idx_w = tmp.idx_w;
if(new_curves==1)
  w= load(sprintf('%s/%s/data/%s/%s/%s/w.dat', directory_name, Fig_num,experiment_name,region, bound_type));
  w_rand = w(idx_w,:); 
  save(sprintf('%s/%s/data/%s/%s/%s/w_rand.mat', directory_name, Fig_num,experiment_name, region, bound_type), 'w_rand');
else
  tmp =  load(sprintf('%s/%s/data/%s/%s/%s/w_rand.mat', directory_name, Fig_num,experiment_name, region, bound_type));
  w_rand = tmp.w_rand;
    
end  
w_BACK = load(sprintf('%s/%s/data/%s/%s/%s/w_BACK.dat',  directory_name,Fig_num, experiment_name, region, bound_type));

%%
Vbin = load(sprintf('%s/%s/data/%s/%s/%s/Vbin.dat',directory_name, Fig_num, experiment_name, region, bound_type));
Vspk = load(sprintf('%s/%s/data/%s/%s/%s/Vspk.dat',directory_name, Fig_num, experiment_name, region, bound_type));


SPB_tonic= load(sprintf('%s/%s/data/%s/%s/%s/SPB_tonic.dat',directory_name, Fig_num, experiment_name, region, bound_type));
SPB_burst= load(sprintf('%s/%s/data/%s/%s/%s/SPB_burst.dat',directory_name, Fig_num, experiment_name, region, bound_type));
PER_tonic= load(sprintf('%s/%s/data/%s/%s/%s/PER_tonic.dat',directory_name, Fig_num, experiment_name, region, bound_type));
PER_burst= load(sprintf('%s/%s/data/%s/%s/%s/PER_burst.dat',directory_name, Fig_num, experiment_name, region, bound_type));
IBF_tonic= load(sprintf('%s/%s/data/%s/%s/%s/IBF_tonic.dat',directory_name, Fig_num, experiment_name, region, bound_type));
IBF_burst= load(sprintf('%s/%s/data/%s/%s/%s/IBF_burst.dat',directory_name, Fig_num, experiment_name, region, bound_type));
DC_tonic= load(sprintf('%s/%s/data/%s/%s/%s/DC_tonic.dat',directory_name, Fig_num, experiment_name, region, bound_type));
DC_burst= load(sprintf('%s/%s/data/%s/%s/%s/DC_burst.dat',directory_name, Fig_num, experiment_name, region, bound_type));
freq_tonic= load(sprintf('%s/%s/data/%s/%s/%s/freq_tonic.dat',directory_name, Fig_num, experiment_name, region, bound_type));
freq_burst= load(sprintf('%s/%s/data/%s/%s/%s/freq_burst.dat',directory_name, Fig_num, experiment_name, region, bound_type));








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
color_dep = [246 157 135]./255; 
color_pot = [209 62 69]./255;

%%  PLOT W
points = [1 19 32 47 80] ;

for idx_state=1:1:size(w_BACK,2)
    figure
    hold on
    plot(w_rand(:,(idx_state-1)*Tstate+1:idx_state*Tstate)', 'color', [0.75 0.75 0.75], 'linewidth', 0.25)
    for idx_p=1:1:length(points)
        plot( w_rand(points(idx_p),(idx_state-1)*Tstate+1:idx_state*Tstate)','color', color_p(idx_p, :),'linewidth', 2)
    end
    switch idx_state
        case 1
            state="tonic";
            T1=0;
            T2=25000; 
            T3=50000; 
        case 2
            state="burst";
            T1=50000;
            T2=75000; 
            T3=100000; 
    end

    xticks([T1 T2 T3])
    xticklabels({'','','','',''})
    yticks([0 0.5 1])
    yticklabels({'','', ''})


    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 4.5 3.5]);
    print(sprintf('/Users/kathleen/Documents/PhD/2023-reset/%s/fig/%s_%s_%s_w%s',Fig_num, experiment_name, region, bound_type, state), '-depsc', '-painters')

end
%% ZOOM at the beginning of the STATE

for idx_state=1:1:size(w_BACK,2)
    figure
    hold on
    plot(w_rand(:,(idx_state-1)*Tstate+1:(idx_state-1)*Tstate+10000)', 'color', [0.75 0.75 0.75], 'linewidth', 0.25)
    for idx_p=1:1:length(points)
        plot( w_rand(points(idx_p),(idx_state-1)*Tstate+1:(idx_state-1)*Tstate+10000)','color', color_p(idx_p, :),'linewidth', 2)
    end
    switch idx_state
        case 1
            state="tonic";
            T1=0;
            T2=5000; 
            T3=10000; 
        case 2
            state="burst";
            T1=50000;
            T2=55000; 
            T3=60000; 
    end

    xticks([T1 T2 T3])
    xticklabels({'','','','',''})
    yticks([0 0.5 1])
    yticklabels({'','', ''})


    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 5 3.5]);
    print(sprintf('/Users/kathleen/Documents/PhD/2023-reset/%s/fig/%s_%s_%s_w_ZOOMstart%s',Fig_num, experiment_name, region, bound_type, state), '-depsc', '-painters')

end
%% ZOOM At the end of the simulation


for idx_state=1:1:size(w_BACK,2)
    figure
    hold on
    plot(w_rand(:,idx_state*Tstate-10000+1:idx_state*Tstate)', 'color', [0.75 0.75 0.75], 'linewidth', 0.25)
    for idx_p=1:1:length(points)
        plot( w_rand(points(idx_p),idx_state*Tstate-10000+1:idx_state*Tstate)','color', color_p(idx_p, :),'linewidth', 2)
    end
    switch idx_state
        case 1
            state="tonic";
            T1=40001;
            T2=45000; 
            T3=50000; 
        case 2
            state="burst";
            T1=90001;
            T2=95000; 
            T3=100000; 
    end

    xticks([T1 T2 T3])
    xticklabels({'','','','',''})
    yticks([0 0.5 1])
    yticklabels({'','', ''})


    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 3 3.5]);
    print(sprintf('/Users/kathleen/Documents/PhD/2023-reset/%s/fig/%s_%s_%s_w_ZOOMend%s',Fig_num, experiment_name, region, bound_type, state), '-depsc', '-painters')

end



%%



idx_Vpre = 10; 
idx_Vpost = 53; 
idx_Ca = 19; 

%idx_Vpre = 5; 
%idx_Vpost = 52; 
%idx_Ca = 4; 

for idx_state =1:1:2
    switch idx_state
        case 1
            T1 = Tcycle/2-2500; 
            T2 = Tcycle/2-1500;
        case 2
            T1 = Tcycle/2+1500; 
            T2 = Tcycle/2+2500;
    end
    

    figure
    count=1;
    hold on
    idx_t = find(Vbin(idx_Vpre,:)==1); 
    plot(idx_t, Vbin(idx_Vpre,idx_t)+count,'o','color', 'k', 'MarkerSize', 2)
    count=count-0.01;
    idx_t = find(Vbin(idx_Vpost,:)==1); 
    plot(idx_t, Vbin(idx_Vpre,idx_t)+count,'o','color', 'k', 'MarkerSize', 2)
     


    xlim([T1 T2])
    %ylim([1.95 2 ])
    xticks([0 1000])
    xticklabels({'',''})
    yticks([0.99 1])
    yticklabels({'',''})

    axis off
    box off


    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 8 0.5]);
    print(sprintf('/Users/kathleen/Documents/PhD/2023-reset/%s/fig/%s_%s_%s_Raster_%d',Fig_num, experiment_name, region, bound_type, idx_state), '-depsc', '-painters')



end


%% prediction vs simulation
wpredict = zeros(size(w_BACK)); 

dt=1;

A_p = 0.0096;
A_m = 0.0053; 
tau_p = 16.8;
tau_m = 33.7; 


for idx=1:1:Ncycles*2
    count=1; 


    T1 = (idx-1)*Tcycle/2 +1;
    T2 = idx*Tcycle/2;

    for idx_post =1:1:50
        for idx_pre=1:1:50
            Vpre = Vbin(idx_pre+1,T1:T2); 
            Vpost = Vbin(51+idx_post,T1:T2); 
           [c, lags] = xcorr(Vpost,Vpre);

            A = zeros (2, length(lags));
            A(1,:) = lags(:); % index
            A(2,:) = c(:); % valeurs de la correlation 

            F1=0; 
            F2=0; 

            % calcul de C+
            for s = floor(length(A)/2)+1:1:length(A)

                F1 = F1 + exp(-((A(1,s))*dt)/tau_p)*A(2,s);

            end 
            % calcul de C-
            for s = 1:1:floor(length(A)/2)+1

                F2 = F2 + exp(((A(1,s))*dt)/tau_m)*A(2,s);

            end 
            
            
            Cp(count,idx) = F1; 
            Cm(count,idx) = F2; 

            %if(bound_type=="SB")
                lambda = (A_p/A_m)*(F1/F2);
                wpredict(count,idx)=lambda / (1+lambda);
            %else
                wSAT(count,idx) = (A_p*F1 - A_m*F2)/((T2-T1)*dt);
            %end
            count=count+1; 

        end
    end
end

%%
figure
stem(lags*dt,c./(T2*dt-T1*dt))

%%

for idx=1:1:Ncycles*2
    figure(idx)
    hold on 
    
    %if(bound_type=="SB")
        plot([0 1], [0 1], '-', 'color', [0.5 0.5 0.5])
        plot(wpredict(:,idx), w_BACK(:,idx),'o', 'MarkerSize',1)
         xlim([0 1])
         ylim([0 1])
         yticks([0 1])
         yticklabels({'',''})
         xticks([0 1])
         xticklabels({'',''})
    %else
    %    if(idx==1)
    %        wdiff(:,1)=(w_BACK(:,1)-0.5)./(Tcycle*0.5);
    %    else
    %        wdiff(:,2)=(w_BACK(:,2)-w_BACK(:,1))./(Tcycle*0.5);
    %        
    %    end
    %    
    %    plot(wSAT(:,idx), wdiff(:,idx),'o', 'MarkerSize', 1)
    %end

    
   
    box off
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 ptx_hist pty_hist]);
    print(sprintf('/Users/kathleen/Documents/PhD/2023-reset/%s/fig/%s_%s_%s_PREDICT%d',Fig_num, experiment_name, region, bound_type, idx), '-depsc', '-painters')

   
end


   
 %% histogramme Corr


for idx=1:1:2
    figure(idx)
    histogram(Cp(:,idx),'FaceColor', color_hist)
    %yticks([0 1500])
    %yticklabels({'',''})
    %ylim([0 1500])
    %xlim([-0.1 1.1])
    box off
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 ptx_hist pty_hist]);
    print(sprintf('/Users/kathleen/Documents/PhD/2023-reset/%s/fig/%s_%s_%s_Cp%d',Fig_num, experiment_name, region, bound_type, idx), '-depsc', '-painters')
end
%%
for idx=1:1:2
    figure(idx)
    histogram(Cm(:,idx),'FaceColor', color_hist)
    %yticks([0 1500])
    %yticklabels({'',''})
    %ylim([0 1500])
    %xlim([-0.1 1.1])
    box off
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 ptx_hist pty_hist]);
    print(sprintf('/Users/kathleen/Documents/PhD/2023-reset/%s/fig/%s_%s_%s_Cm%d',Fig_num, experiment_name, region, bound_type, idx), '-depsc', '-painters')
end

%% compute square error

SQE = zeros(1,2); 
SQE(1) = immse(wpredict(:,1),w_BACK(:,1)); 
SQE(2) = immse(wpredict(:,2), w_BACK(:,2)); 

csvwrite(sprintf('/Users/kathleen/Documents/PhD/2023-reset/%s/fig/%s_%s_%s_SQE',Fig_num, experiment_name, region, bound_type), SQE)



