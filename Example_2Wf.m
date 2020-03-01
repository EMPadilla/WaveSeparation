
load('DATA_2Wf')

x = Data.x;
d = Data.d;
eta = Data.eta;
etaNoise = Data.etaNoise;
time = Data.time;
fs = Data.fs;
f1 = Data.f1;
Th_IFW = Data.IFW;
Th_OFW = Data.OFW;

%% Cross-shore amplitudes
[Amp.Total,~] = FUNCTION_AmplitudePhase(eta,time,f1); close
[Amp.Th_IFW,~] = FUNCTION_AmplitudePhase(Th_IFW,time,f1); close
[Amp.Th_OFW,~] = FUNCTION_AmplitudePhase(Th_OFW,time,f1); close

%% Propagation characteristics:

% Propagation characteristics for the free waves at f1
[Kinematic.Free.f1] = FUNCTION_Kinematic_FreeWaves(x,d,f1);

%% Quantitative Analysis: Wave-Separation    
Vbound = Kinematic.Free.f1.c*2;
[IBW.P4_Psep3,OBW.P4_Psep3,OFW.P4_Psep3,IFW.P4_Psep3] = FUNCTION_WaveSeparation4W(etaNoise,x,d,Vbound,f1,fs,[4 3]);
[IBW.P4_Psep6,OBW.P4_Psep6,OFW.P4_Psep6,IFW.P4_Psep6] = FUNCTION_WaveSeparation4W(etaNoise,x,d,Vbound,f1,fs,[4 6]);

%% Theoretical estimation of the nodes and antinodes 
WaveNumbers(1,:) = Kinematic.Free.f1.Keq * cos(0);
WaveNumbers(2,:) = Kinematic.Free.f1.Keq * cos(pi); 
[xAnti.IFW_OFW,xNodes.IFW_OFW,L_Approx.IFW_OFW] = FUNCTION_xAntinodes2D(x,[IFW.P4_Psep6.Ph0 OFW.P4_Psep6.Ph0],WaveNumbers);


%% Latex Interpreter
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaulttextinterpreter','latex');
set(0, 'defaultColorbarTickLabelInterpreter','latex');

%% Graph Qualitative Analysis

EtaColor  = [0 0.447058826684952 0.74117648601532;0.03125 0.46433824300766 0.749264717102051;0.0625 0.481617659330368 0.757352948188782;0.09375 0.498897075653076 0.765441179275513;0.125 0.516176462173462 0.773529410362244;0.15625 0.533455908298492 0.781617641448975;0.1875 0.550735294818878 0.789705872535706;0.21875 0.568014681339264 0.797794103622437;0.25 0.585294127464294 0.805882334709167;0.28125 0.60257351398468 0.813970625400543;0.3125 0.619852960109711 0.822058856487274;0.34375 0.637132346630096 0.830147087574005;0.375 0.654411792755127 0.838235318660736;0.40625 0.671691179275513 0.846323549747467;0.4375 0.688970565795898 0.854411780834198;0.46875 0.706250011920929 0.862500011920929;0.5 0.723529398441315 0.87058824300766;0.53125 0.740808844566345 0.878676474094391;0.5625 0.758088231086731 0.886764705181122;0.59375 0.775367677211761 0.894852936267853;0.625 0.792647063732147 0.902941167354584;0.65625 0.809926450252533 0.911029398441315;0.6875 0.827205896377563 0.919117629528046;0.71875 0.844485282897949 0.927205860614777;0.75 0.86176472902298 0.935294151306152;0.78125 0.879044115543365 0.943382382392883;0.8125 0.896323502063751 0.951470613479614;0.84375 0.913602948188782 0.959558844566345;0.875 0.930882334709167 0.967647075653076;0.90625 0.948161780834198 0.975735306739807;0.9375 0.965441167354584 0.983823537826538;0.96875 0.982720613479614 0.991911768913269;1 1 1;0.99519294500351 0.978241622447968 0.970904469490051;0.990385830402374 0.956483244895935 0.941808998584747;0.985578775405884 0.934724867343903 0.912713468074799;0.980771660804749 0.91296648979187 0.88361793756485;0.975964605808258 0.891208112239838 0.854522466659546;0.971157491207123 0.869449734687805 0.825426936149597;0.966350436210632 0.847691357135773 0.796331465244293;0.961543321609497 0.82593297958374 0.767235934734344;0.956736266613007 0.804174602031708 0.738140404224396;0.951929152011871 0.782416224479675 0.709044933319092;0.947122097015381 0.760657787322998 0.679949402809143;0.942314982414246 0.738899409770966 0.650853872299194;0.937507927417755 0.717141032218933 0.62175840139389;0.93270081281662 0.695382654666901 0.592662870883942;0.927893757820129 0.673624277114868 0.563567340373993;0.923086643218994 0.651865899562836 0.534471869468689;0.918279588222504 0.630107522010803 0.50537633895874;0.913472473621368 0.608349144458771 0.476280838251114;0.908665418624878 0.586590766906738 0.447185337543488;0.903858304023743 0.564832389354706 0.418089807033539;0.899051249027252 0.543074011802673 0.388994306325912;0.894244134426117 0.521315634250641 0.359898805618286;0.889437079429626 0.499557256698608 0.33080330491066;0.884629964828491 0.477798879146576 0.301707774400711;0.879822909832001 0.456040501594543 0.272612273693085;0.875015795230865 0.434282094240189 0.243516758084297;0.870208740234375 0.412523716688156 0.214421257376671;0.86540162563324 0.390765339136124 0.185325741767883;0.860594570636749 0.369006961584091 0.156230241060257;0.855787456035614 0.347248584032059 0.127134725451469;0.850980401039124 0.325490206480026 0.0980392172932625];
Fig1 = figure('Colormap',EtaColor,'units','points','position',[50,50,950,400]);

ax1 = subplot(2,2,1); 
    clear p
    p(5) = area([xAnti.IFW_OFW(1); xAnti.IFW_OFW(1)+L_Approx.IFW_OFW],[0.02; 0.02],0,'LineStyle','none');hold on; p(5).FaceColor = [0.85 0.85 0.85];hold on
    p(1) = plot(x,Amp.Total,'.:','Color',[0 0.4470 0.7410],'markersize',8);hold on
    p(2) =plot(x,Amp.Th_OFW,'r-');
    p(3) =plot(x,Amp.Th_IFW,'g-');
    p(4) = plot(xAnti.IFW_OFW,interp1(x,Amp.Total,xAnti.IFW_OFW),'r.','markersize',22);
    legend(p,'Total','OFW','IFW','$X_{anti}[IFW,OFW]$','$L_{IFW,OFW}$','Location','southeastoutside')
    ylim([0 0.012])
    xlim([x(1) x(end)])
    ylabel('Wave Amplitude (m)')
    set(gca,'fontsize',12)
    grid on 
    text(0.5,0.011,'a.','FontSize',18)

ax2 = subplot(2,2,2);
    clear p
    [meshX,meshT] = meshgrid(x,time);
    contour(meshX(1:1:1000,:),meshT(1:1:1000,:),eta(1:1000,:),'LevelStep',Amp.Total(end)/10,'Fill','on');hold on;
    c = colorbar('Location','eastoutside');
    ylabel(c,'$\eta_{2f_1}$ (m)','Interpreter','latex');
    p(1) = plot(x,Kinematic.Free.f1.tau-1,'g-','LineWidth',2);hold on
    p(2) = plot(x,-Kinematic.Free.f1.tau+11.6,'r-','LineWidth',2);hold on
    legend(p,'IFW','OFW','Location','West')
    ylim([0 12])
    xlim([x(1) x(end)])
    ylabel('Time (s)')
    set(gca,'fontsize',12)
    grid on
    text(0.5,10.7,'c.','FontSize',18)    
    
ax3 = subplot(2,2,3);    
    clear p
    p(2) = plot(x,eta(163,:),'k-');hold on
    p(3) = plot(x,eta(200,:),'k--');hold on
    p(1) = plot(x,Amp.Total,'-','LineWidth',2);hold on
    legend(p,'$\hat{\eta}$',sprintf('t = %ss',num2str(time(163))),sprintf('t = %ss',num2str(time(200))),'Location','southeastoutside')
    ylabel('$\eta$ (m)')
    ylim([-0.012 0.012])
    xlim([x(1) x(end)])
    xlabel('$x$ (m)')
    set(gca,'fontsize',12)
    grid on
    text(0.5,0.009,'b.','FontSize',18)

ax4 = subplot(2,2,4);
    plot(x,-d,'k-')
    ylim([-0.6 0])
    xlim([x(1) x(end)])
    xlabel('$x$ (m)')
    ylabel('Depth (m)')
    set(gca,'fontsize',12)
    grid on
    text(0.5,-0.2,'d.','FontSize',18)  
    
    
set(ax1,'units','points','position',[50,50+100+50,300,150]) 
set(ax2,'units','points','position',[50+500,50+50+50,300,200])
set(ax3,'units','points','position',[50,50,300,100]) 
set(ax4,'units','points','position',[50+500,50,300,50])   
clear ax1 ax2 ax3 ax4

    
%% Graph Quantitative Analysis: Wave-Separation    

Fig2 = figure('units','points','position',[0,0,800,450]);
    
ax1 = subplot(2,2,1);
    clear p
    p(1) = plot(x,Amp.Total,'.:','Color',[0 0.4470 0.7410],'markersize',8);hold on
    p(2) = plot(x,IBW.P4_Psep3.Amp,'.k','markersize',10);
    p(3) = plot(x,OBW.P4_Psep3.Amp,'.','Color',[0.93 0.69 0.13],'markersize',8); 
    p(4) = plot(x,OFW.P4_Psep3.Amp,'.r','markersize',8);
    p(5) = plot(x,Amp.Th_OFW,'r--');
    p(6) = plot(x,IFW.P4_Psep3.Amp,'.g','markersize',8);
    p(7) = plot(x,Amp.Th_IFW,'--g');
    legend(p,'Total','IBW (Sep)','OBW (Sep)','OFW (Sep)','OFW (Th)','IFW (Sep)','IFW (Th)','Location','north')
    xlim([x(1) x(end)])
    ylim([0 0.012])
    ylabel('Wave Amplitude (m)')
    grid on
    set(gca,'fontsize',12)
    text(0.75,0.011,'a.','FontSize',18)
    
ax2 = subplot(2,2,2);
    plot(IBW.P4_Psep3.alpha,'.:','Color',[0 0.4470 0.7410],'markersize',12);hold on
    ylabel('$\alpha$')
    xlim([1 6])
    grid on
    set(gca,'fontsize',12)
    text(1.3,0.7,'b.','FontSize',18)
    
ax3 = subplot(2,2,3);
    clear p
    p(1) = plot(x,Amp.Total,'.:','Color',[0 0.4470 0.7410],'markersize',8);hold on
    p(2) = plot(x,IBW.P4_Psep6.Amp,'.k','markersize',10);
    p(3) = plot(x,OBW.P4_Psep6.Amp,'.','Color',[0.93 0.69 0.13],'markersize',8);
    p(4) = plot(x,OFW.P4_Psep6.Amp,'.r','markersize',8);
    p(5) = plot(x,Amp.Th_OFW,'r--');
    p(6) = plot(x,IFW.P4_Psep6.Amp,'.g','markersize',8);
    p(7) = plot(x,Amp.Th_IFW,'--g');
    legend(p,'Total','IBW (Sep)','OBW (Sep)','OFW (Sep)','OFW (Th)','IFW (Sep)','IFW (Th)','Location','north')
    xlim([x(1) x(end)])
    ylim([0 0.012])
    ylabel('Wave Amplitude (m)')
    grid on
    set(gca,'fontsize',12)
    text(0.75,0.011,'c.','FontSize',18)

ax4 = subplot(2,2,4);
    plot(IBW.P4_Psep6.alpha,'.:','Color',[0 0.4470 0.7410],'markersize',12);hold on
    ylabel('$\alpha$')
    xlim([1 6])
    grid on
    set(gca,'fontsize',12)
    text(1.3,0.7,'d.','FontSize',18)    
    
set(ax1,'units','points','position',[50,50+150+50,350,150]) 
set(ax2,'units','points','position',[50+400,50+150+50,150,150])
set(ax3,'units','points','position',[50,50,350,150]) 
set(ax4,'units','points','position',[50+400,50,150,150])   
clear ax1 ax2 ax3 ax4    
    