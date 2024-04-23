clear all
close all
clc

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

directory_name = "/Users/kathleen/Documents/PhD/2023-Project"; 
experiment_name = "Graupner2016"; 
bound_type = "HB";
region = "Cortex"; 
Fig_num= "Fig3_Demo"; 
new_curves=1; 

Ncycles=1; 
Tcycle=100000; 
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
Ca = load(sprintf('%s/%s/data/%s/%s/%s/Ca.dat', directory_name,Fig_num, experiment_name, region, bound_type));


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

%%
demo_pot = load(sprintf('%s/%s/data/%s/%s/%s/demo_pot.dat',directory_name, Fig_num, experiment_name, region, bound_type));
demo_dep = load(sprintf('%s/%s/data/%s/%s/%s/demo_dep.dat',directory_name, Fig_num, experiment_name, region, bound_type));
demo_nul = load(sprintf('%s/%s/data/%s/%s/%s/demo_nul.dat',directory_name, Fig_num, experiment_name, region, bound_type));


demo_pot = demo_pot./(Tcycle*0.5/dt);
demo_dep = demo_dep./(Tcycle*0.5/dt);
demo_nul = demo_nul./(Tcycle*0.5/dt);

%%

    wMAX=1;
    wMIN=0; 
    Tstate=50000; 
    
    demo_lim= zeros(size(w,1),Ncycles*2);
    for idx_state = 1:2:Ncycles*2
        for idw = 1:1:size(w,1)
            if(w(idw,idx_state*Tstate)>=wMAX)
                idx_win = find(w(idw, (idx_state-1)*Tstate+100:idx_state*Tstate)>=wMAX); 
                demo_lim(idw, idx_state)=(idx_state-1)*Tstate+100+idx_win(1);
            else
                if (w(idw,idx_state*Tstate)<=wMIN)
                    idx_win = find(w(idw, (idx_state-1)*Tstate+100:idx_state*Tstate)<=wMIN); 
                    demo_lim(idw, idx_state)=(idx_state-1)*Tstate+100+idx_win(1);
                else
                    demo_lim(idw,idx_state)=idx_state*Tstate; 
                end
            end
            if(w(idw,idx_state*Tstate)==w(idw, (idx_state-1)*Tstate+100))
                demo_lim(idw,idx_state)=0; 
            end

        end
    end


    for idx_state = 2:2:Ncycles*2
        for idw = 1:1:size(w,1)
            if(w(idw,idx_state*Tstate)>=wMAX)
                idx_win = find(w(idw, (idx_state-1)*Tstate+1000:idx_state*Tstate)>=wMAX); 
                demo_lim(idw, idx_state)=(idx_state-1)*Tstate+1000+idx_win(1);
            else
                if (w(idw,idx_state*Tstate)<=wMIN)
                    idx_win = find(w(idw, (idx_state-1)*Tstate+1000:idx_state*Tstate)<=wMIN); 
                    demo_lim(idw, idx_state)=(idx_state-1)*Tstate+1000+idx_win(1);
                else
                    demo_lim(idw,idx_state)=idx_state*Tstate; 
                end
            end 
            if(w(idw,idx_state*Tstate)==w(idw, (idx_state-1)*Tstate+500))
                demo_lim(idw,idx_state)=0; 
            end

        end
    end
    
    
%%


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
color_dep = [255 192 174]./255;%[246 157 135]./255; 
color_pot = [209 62 69]./255; %[146 196 196]./255; %

%% w
points = [1 19 32 47 80] ;
figure
hold on
plot(w_rand', 'color', [0.75 0.75 0.75], 'linewidth', 0.25)
for idx_p=1:1:length(points)
    plot( w_rand(points(idx_p),:)','color', color_p(idx_p, :),'linewidth', 2)
end

xticks([0 25000 50000 75000 100000])
xticklabels({'','','','',''})
yticks([0 0.5 1])
yticklabels({'','', ''})


set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 ptx_long pty]);
print(sprintf('/Users/kathleen/Documents/PhD/2023-Project/%s/fig/%s_%s_%s_w',Fig_num, experiment_name, region, bound_type), '-depsc', '-painters')


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
            T1=1;
            T2=5000; 
            T3=10000; 
        case 2
            state="burst";
            T1=1;
            T2=5000; 
            T3=10000; 
    end

    xticks([T1 T2 T3])
    xticklabels({'','','','',''})
    yticks([0 0.5 1])
    yticklabels({'','', ''})


    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 5 3.5]);
    print(sprintf('/Users/kathleen/Documents/PhD/2023-Project/%s/fig/%s_%s_%s_w_ZOOMstart%s',Fig_num, experiment_name, region, bound_type, state), '-depsc', '-painters')

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
            T1=1;
            T2=5000; 
            T3=10000; 
        case 2
            state="burst";
            T1=1;
            T2=5000; 
            T3=10000; 
    end

    xticks([T1 T2 T3])
    xticklabels({'','','','',''})
    yticks([0 0.5 1])
    yticklabels({'','', ''})


    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 3 3.5]);
    print(sprintf('/Users/kathleen/Documents/PhD/2023-Project/%s/fig/%s_%s_%s_w_ZOOMend%s',Fig_num, experiment_name, region, bound_type, state), '-depsc', '-painters')

end



%% SLOPE NUMERIQUE

slope = zeros(size(w_BACK,1),Ncycles*2);
nPre= 50; 
nPost=50; 
Tstate =Tcycle/2; 
w0 = w(:,2);
for idx=1:1:Ncycles*2
    for idw = 1:1:size(slope,1)
        slope(idw, idx) = (w_BACK(idw, idx) - w0(idw))./ (demo_lim(idw, idx) - (idx-1)*Tstate); 
        if(demo_lim(idw,idx)==0)
            slope(idw,idx) =0; 
        end
    end
    if(idx<Ncycles*2)
        w0 = w(:, idx*Tstate+1);
    end
end

%% reconstruction

if(new_curves==1)
    figure
    hold on
    plot(w_rand', 'color', [0.75 0.75 0.75], 'linewidth', 0.25)


    for idx_p=1:1:length(points)
        plot( w_rand(points(idx_p),:)','color', color_p(idx_p, :),'linewidth', 2)
        w0 = w(:,2);
        for idx_state=1:1:Ncycles*2
            plot( [(idx_state-1)*Tstate, demo_lim(idx_w(points(idx_p)),idx_state)], [w0(idx_w(points(idx_p))) w0(idx_w(points(idx_p)))+slope(idx_w(points(idx_p)),idx_state).*(demo_lim(idx_w(points(idx_p)),idx_state)-(idx_state-1)*Tstate)], '-o','color', color_p(idx_p, :),'linewidth', 2)
            if(idx_state<Ncycles*2)
                w0 = w(:,idx_state*Tstate+1);
            end
        end

    end
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
    print(sprintf('/Users/kathleen/Documents/PhD/2023-Project/%s/fig/%s_%s_%s_Raster_%d',Fig_num, experiment_name, region, bound_type, idx_state), '-depsc', '-painters')



    % Hypnogramm

    %T1 = Tcycle/2-2000; 
    %T2 = Tcycle/2-1000;

    colors = zeros(3,3); 
    colors(1,:) = [ 1 1 1]; 
    colors(2,:) = color_dep; 
    colors(3,:) = color_pot; 

    thresholds = [1 2]; 
     
    figure
    hold on
    plot(T1:1:T2, Ca(idx_Ca,T1:T2), 'color', [0.3 0.3 0.3]) 
    plot([T1 T2], [1 1], 'color', color_dep)
    plot([T1 T2], [2 2], 'color', color_pot)
    xlim([T1 T2])
    ylim([0 4])
    xticks([0 1000])
    xticklabels({'',''})
    %yticks([0.95 1])
    yticklabels({'',''})
    axis off
    box off
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 8 1.5]);
    print(sprintf('/Users/kathleen/Documents/PhD/2023-Project/%s/fig/%s_%s_%s_Ca%d_%d',Fig_num,experiment_name, region, bound_type, idx_Ca, idx_state), '-depsc', '-painters')
    
    

    figure
    set_background(T1:1:T2,Ca(idx_Ca,T1:T2),colors,thresholds)
    xlim([T1 T2])
    ylim([0 1])
    xticks([0 1000])
    xticklabels({'',''})
    %yticks([0.95 1])
    yticklabels({'',''})
    axis off
    box off
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 8 1.5]);
    print(sprintf('/Users/kathleen/Documents/PhD/2023-Project/%s/fig/%s_%s_%s_Hypnogram%d_%d',Fig_num,experiment_name, region, bound_type, idx_Ca , idx_state), '-depsc', '-painters')

    figure
    T1 = 0;
    T2 = demo_pot(idx_Ca,idx_state); 
    T3 = demo_pot(idx_Ca,idx_state)+demo_dep(idx_Ca,idx_state); 
    T4 = demo_pot(idx_Ca,idx_state)+demo_dep(idx_Ca,idx_state)+demo_nul(idx_Ca, idx_state); 
    vPot  = [T1 0;T2  0; T2 1; T1 1];
    vDep = [T2 0; T3 0;T3 1; T2 1];
    vNul = [T3 0; T4 0; T4 1; T3 1]; 
    f = [1 2 3 4];
    patch('Faces',f,'Vertices', vPot, 'FaceColor', color_pot, 'EdgeColor', 'none','FaceAlpha', 1)
    patch('Faces',f,'Vertices', vDep, 'FaceColor', color_dep, 'EdgeColor', 'none','FaceAlpha', 1)
    patch('Faces',f,'Vertices', vNul, 'FaceColor', [1 1 1], 'EdgeColor', 'none','FaceAlpha', 1)
    
    xlim([0 1])
    ylim([0 1])
    xticks([0 1])
    xticklabels({'',''})
    yticks([0.95 1])
    yticklabels({'',''})
    axis off
    box off
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 8 1.5]);
    print(sprintf('/Users/kathleen/Documents/PhD/2023-Project/%s/fig/%s_%s_%s_Cumuls%d_%d',Fig_num,experiment_name, region, bound_type, idx_Ca , idx_state), '-depsc', '-painters')

    
end



%% prediction vs simulation

tau_w    = 520.76129e3;
gamma_p  = 597.08922;
gamma_d  = 137.7586;
theta_p = 2.009289;
theta_d = 1.0;

if(bound_type=="SB")
    Omega_p = gamma_p/(gamma_p+gamma_d);
    Omega_d = 0;
    Omega_0 = 0;

    tauw_p = tau_w /(gamma_p + gamma_d);
    tauw_d = tau_w / gamma_d;
    tauw_0 = 0;
    zeta = tauw_d / tauw_p;
else
    Omega_p =gamma_p*0.5-gamma_d*0.5;
    Omega_d = -0.5*gamma_d;
    Omega_0 = 0;

    tauw_p = tau_w;%/(gamma_p + gamma_d);
    tauw_d = tau_w;%/ gamma_d;
    tauw_0 = 0;
    zeta = tauw_d / tauw_p;
DT=3000;
end

slope_predict = zeros(size(demo_pot)); 
count=1;
Mks=4.5; 
%%

%wpredict(:,1) = (Omega_p*demo_pot(:,1)*1/tauw_p+Omega_d*demo_dep(:,1)*1/tauw_d)./(demo_pot(:,1)*1/tauw_p+demo_dep(:,1)*1/tauw_d); 
%wpredict(:,2) = (Omega_p*demo_pot(:,2)*1/tauw_p+Omega_d*demo_dep(:,1)*1/tauw_d)./(demo_pot(:,2)*1/tauw_p+demo_dep(:,2)*1/tauw_d); 

slope_predict(:,1) = (demo_pot(:,1)*1/tauw_p*Omega_p)+(Omega_d*demo_dep(:,1)*1/tauw_d);
slope_predict(:,2) = (demo_pot(:,2)*1/tauw_p*Omega_p)+(Omega_d*demo_dep(:,2)*1/tauw_d);

%slope_predict(:,1) = (demo_pot(:,1)*0.5*(gamma_p-gamma_d)*1/tau_w)-0.5*demo_dep(:,1)*gamma_d*1/tau_w;
%slope_predict(:,2) = (demo_pot(:,2)*1/tauw_p*Omega_p)+(Omega_d*demo_dep(:,2)*1/tauw_d);

for idw=1:1:size(w,1)
    for idx_state=1:1:2
        if(round(demo_nul(idw,idx_state),4)==1)
            slope_predict(idw,idx_state)=0;%w(idw,(idx_state-1)*Tstate+1);
        end
    end
    
end

%% compute square error

SQE = zeros(1,2); 
SQE(1) = immse(slope_predict(:,1),slope(:,1))
SQE(2) = immse(slope_predict(:,2), slope(:,2))


csvwrite(sprintf('/Users/kathleen/Documents/PhD/2023-Project/%s/fig/%s_%s_%s_SQE',Fig_num, experiment_name, region, bound_type), SQE)



%%


for idx=1:1:Ncycles*2
    figure(idx)
    hold on 
    
        %plot([0 1], [0 1], '-', 'color', [0.5 0.5 0.5])
        plot(slope_predict(:,idx), w_BACK(:,idx),'o','Markersize', 1,'MarkerEdgeColor',color_hist,'MarkerFaceColor',color_hist); 
        for idx_p=1:1:length(points)
            plot(slope_predict(idx_w(points(idx_p)),idx)', slope(idx_w(points(idx_p),idx))','o', 'Markersize', Mks,'MarkerEdgeColor',color_p(idx_p,:),'MarkerFaceColor',color_p(idx_p,:)); 
        end
        %xlim([0 1])
        % ylim([0 1])
        % yticks([0 1])
        % yticklabels({'',''})
        % xticks([0 1])
        % xticklabels({'',''})


    
   
    box off
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 ptx_hist pty_hist]);
    print(sprintf('/Users/kathleen/Documents/PhD/2023-Project/%s/fig/%s_%s_%s_PREDICT%d',Fig_num, experiment_name, region, bound_type, idx), '-depsc', '-painters')

   
end



%% histogramme alpha

edges = [0:0.05:1];


for idx_state=1:1:2
    switch idx_state
        case 1
            STATE='TONIC';
        case 2
            STATE='BURST';
    end
    for idx=1:1:3
        switch idx
            case 1
                demo_plot = demo_pot(:,idx_state);
            case 2
                demo_plot = demo_dep(:,idx_state);
            case 3
                demo_plot = demo_nul(:,idx_state);
        end
        figure
        histogram(demo_plot,'FaceColor', color_hist, 'BinEdges', edges)
        yticks([0 1000])
        yticklabels({'',''})
        xticks([0 1])
        xticklabels({'',''})
        ylim([0 1200])
        xlim([-0.1 1.1])
        box off
        set(gcf,'PaperPositionMode','auto');
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 ptx_hist pty_hist]);
        print(sprintf('/Users/kathleen/Documents/PhD/2023-Project/%s/fig/%s_%s_%s_%s%d',Fig_num, experiment_name, region, bound_type, STATE,  idx), '-depsc', '-painters')

        
    end
end    

