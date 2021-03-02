%% FIGURE 1: continuity stimuli schematic
clear all
close all

FS = 12;

figure(1), clf

P = [0.15    0.52    0.22    0.3;
     0.45    0.52   0.22    0.3;
     0.75    0.52     0.22    0.3;
     0.15    0.1    0.22    0.3;
     0.45    0.1    0.22    0.3;
     0.75    0.1    0.22    0.3];
 
t = linspace(0,2,700);
f = 10;
tStart = .2; tStop = 1.8; tGapStart = .8; tGapStop = 1.2;
subplot('position',P(1,:))
    plot(t,(t>=tStart).*(t<=tStop).*((t<=tGapStart) | (t>=tGapStop)).*sin(2*pi*f*(t-tStart)),'k','linewidth',1)
    axis off
    ylim([-1.5 1.5])
    text(-1.22,0,'stimulus','fontsize',FS)
    title('A')
    set(gca,'fontsize',FS)
    
subplot('position',P(2,:))
    x = (t>=tStart).*(t<=tStop).*((t<=tGapStart) | (t>=tGapStop)).*sin(2*pi*f*(t-tStart));
    x((t>=tGapStart) & (t<=tGapStop)) = nan;
    y = (t>=tGapStart).*(t<=tGapStop).*(1*(rand(size(t))-.5));
    y(t<=tGapStart) = nan; y(t>=tGapStop) = nan;
    plot(t,x,'k','linewidth',1)
    hold all
    plot(t,y,'color',[1 1 1]*.5,'linewidth',1)
    axis off
    ylim([-1.5 1.5])
    title('B')
    set(gca,'fontsize',FS)
    
subplot('position',P(3,:))
    x = (t>=tStart).*(t<=tStop).*((t<=tGapStart) | (t>=tGapStop)).*sin(2*pi*f*(t-tStart));
    x((t>=tGapStart) & (t<=tGapStop)) = nan;
    y = (t>=tGapStart).*(t<=tGapStop).*(3*(rand(size(t))-.5));
    y(t<=tGapStart) = nan; y(t>=tGapStop) = nan;
    plot(t,x,'k','linewidth',1)
    hold all
    plot(t,y,'color',[1 1 1]*.5,'linewidth',1)
    axis off
    ylim([-1.5 1.5])
    title('C')
    set(gca,'fontsize',FS)

subplot('position',P(4,:))
    x = 0*(t>=tStart).*(t<=tStop).*((t<=tGapStart) | (t>=tGapStop));
    x((t>=tGapStart) & (t<=tGapStop)) = nan;
    x(t<tStart) = nan; x(t>tStop) = nan;
    plot(t,x,'k','linewidth',6)
    text(-1.3,.3,'perceived','fontsize',FS)
    text(-1.05,-.3,'tone','fontsize',FS)
    axis off
    ylim([-1 1])

subplot('position',P(5,:))
    x = 0*(t>=tStart).*(t<=tStop).*((t<=tGapStart) | (t>=tGapStop));
    x((t>=tGapStart) & (t<=tGapStop)) = nan;
    x(t<tStart) = nan; x(t>tStop) = nan;
    plot(t,x,'k','linewidth',6)
    axis off
    ylim([-1 1])

    
subplot('position',P(6,:))
    x = 0*(t>=tStart).*(t<=tStop).*((t<=tGapStart) | (t>=tGapStop));
    x(t<tStart) = nan; x(t>tStop) = nan;
    plot(t,x,'k','linewidth',6)
    axis off
    ylim([-1 1])

set(gcf,'units','centimeters','position',[0 0 18 5])
set(gcf, 'PaperPositionMode','auto') 
 
    
%% FIGURE 2: 6 regions. equilibrium curves 


clear all
close all

FS = 12;

%locations of knees to define three models
IkneeBistable = [-2 2]; % BISTABLE MODEL %
IkneeHysteresis = [.2 1]; % HYSTERESIS MODEL %
IkneeCombined = [.2 6]; % COMBINED MODEL %

%%% aE and m calculator
k = 1;
xR = @(aE) 1/2*(1-sqrt(1-4*k./aE)); % right knee at rest
xL = @(aE) 1/2*(1+sqrt(1-4*k./aE)); % lower knee at rest
aEm = @(I) fsolve(@(P) [I(1)-P(2)+log(1/xL(P(1))-1)+P(1)*xL(P(1)) ,  I(2)-P(2)+log(1/xR(P(1))-1)+P(1)*xR(P(1))  ]  , [5 5]);

%%% aE and m calculations
Iknee = IkneeHysteresis; % hysteresis
a = aEm([Iknee]); AE(1) = a(1); M(1) = a(2);

Iknee = IkneeBistable; % bistable
a = aEm([Iknee]); AE(2) = a(1); M(2) = a(2);

Iknee = IkneeCombined; % combined
a = aEm([Iknee]); AE(3) = a(1); M(3) = a(2);

x=[ linspace(0,.0001,1e3), linspace(.0001,.99,1e3), linspace(.99,1,1e3)];
Ieq = @(x,m,aE) m - k*log(1./x - 1)-aE*x;

C = lines(3);

P = [0.08    0.2    0.36    0.7;
     0.62    0.2    0.36    0.7];

npt = 300;

for i=1:npt
    % region II/III boundary: right knee =0
    j = linspace(-20,0,npt);
    d =  aEm([j(i) 0]); xR23(i) = d(1); yR23(i) = d(2);
    
    % region III/IV boundary: left knee =0 [bistable boundary]  
    j = linspace(0,10,npt);
    d =  aEm([0 j(i)]); xR34(i) = d(1); yR34(i) = d(2);
    
    % region IV/V boundary: right knee = ITmax    
    j = linspace(0,5,npt);
    d =  aEm([j(i) 5]); xR45(i) = d(1); yR45(i) = d(2);
    
    % region V/VI boundary: left knee = ITmax [not excitable]
    j = linspace(5,7.8,npt);
    d =  aEm([5 j(i)]); xR56(i) = d(1); yR56(i) = d(2);
end

figure(2), clf

    % panel A: equilibrium curves
    subplot('position',P(1,:)), hold all

    for i = 1:3
        I = Ieq(x,M(i),AE(i));
        idR = find(x<=xR(AE(i)),1,'last'); x1 =x(1:idR); I1 = I(1:idR);
        idL = find(x>=xL(AE(i)),1,'first'); x3 =x(idL:end); I3 = I(idL:end);
        x2 = x(idR:idL); I2 = I(idR:idL);
        h(3*(i-1)+1) = plot(I1,x1,'linewidth',3,'color',C(i,:));
        h(3*(i-1)+2) = plot(I2,x2,':','linewidth',3,'color',C(i,:));
        h(3*(i-1)+3) = plot(I3,x3,'linewidth',3,'color',C(i,:));
    end
    
    plot([5 5 ], [0 1],'color',.5*[1 1 1],'linewidth',2)
    axis([-2.5 6.5 0 1])
    set(gca,'xtick',[-2:2:6])
    xlabel('tone input (I_{T})','fontsize',FS)
    ylabel('firing rate (x)','fontsize',FS)
    ti=title('A','fontsize',FS); set(ti,'position',get(ti,'position')+[-5.8 -.02 0])
    
    % panel B: regions of aE
    xA = linspace(4,20);
    xM = linspace(0,15);
    [meshA,meshM] = meshgrid(xA,xM);
    for i=1:size(meshA,1)
        for j = 1:size(meshA,2)
            IR(i,j) = Ieq(xR(meshA(i,j)),meshM(i,j),meshA(i,j)); % left knee
            IL(i,j) = Ieq(xL(meshA(i,j)),meshM(i,j),meshA(i,j)); % right knee
        end
    end
    
    subplot('position',P(2,:)), hold all
    plot([4 4 ], [0 15],'k','linewidth',3)
    plot(xR23,yR23,'color','k','linewidth',3);
    plot(xR45,yR45,'color',.5*[1 1 1],'linewidth',3);
    plot(xR34,yR34,'color',0*[1 1 1],'linewidth',3);
    plot(xR56,yR56,'color','k','linewidth',3);
    plot(AE(1),M(1),'o','markerfacecolor',C(1,:),'markersize',10,'markeredgecolor','none')
    plot(AE(2),M(2),'o','markerfacecolor',C(2,:),'markersize',10,'markeredgecolor','none')
    plot(AE(3),M(3),'o','markerfacecolor',C(3,:),'markersize',10,'markeredgecolor','none')
    text(2,9.7,'I','fontsize',FS)
    text(12.5,1.6,'V','fontsize',FS)
    text(12.5,6.,'IV','fontsize',FS)
    text(9.,9.7,'IIIb','fontsize',FS)
    text(5.,6.,'IIIa','fontsize',FS)
    text(5,9.7,'II','fontsize',FS)
    axis([0 13.5 0 10.5])
    xlabel('recurrent excitation (a_E)','fontsize',FS)
    ylabel('half-max of f (m)','fontsize',FS)
    ti=title('B','fontsize',FS); set(ti,'position',get(ti,'position')+[-8.7 -.1 0])
    
    set(gcf,'units','centimeters','position',[0 0 18 7])
    set(gcf, 'PaperPositionMode','auto') 


%% FIGURE 3: Model 1 [hysteresis + sustained inputs]

clear all
close all

FS = 10;

figure(3)

dx1 = .23; dx2= .155;
P = [0.08    0.63   dx2    0.3;
     .08+dx1    0.63   dx2    0.3;
     .08+2*dx1    0.63   dx2   0.3;
     .1+3*dx1    0.63   dx2-.04    0.3;
     .08    0.12   dx2    0.3;
     .08+dx1    0.12   dx2    0.3;
     .08+2*dx1    0.12   dx2    0.3;
     .1+3*dx1    0.12   dx2    0.3];

 % colors [https://davidmathlogic.com/colorblind/#%23D81B60-%231E88E5-%23FFC107-%23004D40]
 C = [0 90 181; 220 50 32]/255;

 
%locations of knees to define three models
IkneeHysteresis = [.2 1]; % HYSTERESIS MODEL %

INmax = 10 ; ITmax = 5;

%%% aE and m calculator
k = 1;
xR = @(aE) 1/2*(1-sqrt(1-4*k./aE)); % right knee at rest
xL = @(aE) 1/2*(1+sqrt(1-4*k./aE)); % lower knee at rest
aEm = @(I) fsolve(@(P) [I(1)-P(2)+log(1/xL(P(1))-1)+P(1)*xL(P(1)) ,  I(2)-P(2)+log(1/xR(P(1))-1)+P(1)*xR(P(1))  ]  , [5 5]);

%%% aE and m calculations
Iknee = IkneeHysteresis; % hysteresis
a = aEm([Iknee]); aE = a(1); m = a(2);

%%% equilibrium curves 
xeq=[ linspace(0,.0001,1e3), linspace(.0001,.9999,1e3), linspace(.9999,1,1e3)];
Ieq = @(x,m,aE) m - k*log(1./x - 1)-aE*x;

%%% equilibrium curves with noise
xRnoise = @(aI,IN) 1/2*(1-sqrt(1-4*k./(aE + aI*IN))); % right knee with noise
xLnoise = @(aI,IN) 1/2*(1+sqrt(1-4*k./(aE + aI*IN))); % lower knee with noise
IeqNoise = @(x,aI,alph,IN) m - k*log(1./x - 1)-aE*x +  aI*IN*(1-x) -alph*IN;
    
%%% parameters for model 1
aI = 1.139; alph = .181; bet = 0; gon = 0; goff = 0; Xinitial = 0;

%%% stimulus levels
IT = [.5 1.5]; 
INmask = 1; INcont = 8;
tStart = 200; % tone onset

subplot('position',P(1,:)); hold all
    I = Ieq(xeq,m,aE);
    idR = find(xeq<=xR(aE),1,'last'); x1 =xeq(1:idR); I1 = I(1:idR);
    idL = find(xeq>=xL(aE),1,'first'); x3 =xeq(idL:end); I3 = I(idL:end);
    x2 = xeq(idR:idL); I2 = I(idR:idL);
    h(1) = plot(I1,x1,'linewidth',3,'color','k'); h0(1) = h(1);
    h(2) = plot(I2,x2,':','linewidth',3,'color','k');  h0(2) = h(2);
    h(3) = plot(I3,x3,'linewidth',3,'color','k');  h0(3) = h(3);
    Ieq1 = fzero(@(x) Ieq(x,m,aE)-IT(1),0.05); plot(IT(1),Ieq1,'o','markeredge','none','markersize',12,'markerfacecolor',C(1,:))
    Ieq2 = fzero(@(x) Ieq(x,m,aE)-IT(2),0.95); plot(IT(2),Ieq2,'o','markeredge','none','markersize',12,'markerfacecolor',C(2,:))
    ylabel('firing rate (x)')
    axis([-.2 2.2 0 1])
    set(gca,'xtick',[0 :1:7],'ytick',0:.5:1,'fontsize',FS)
    ti = title('A1 (tone)'); 
    set(ti,'position',get(ti,'position')+[-.4 0 0]);

   
subplot('position',P(5,:)), hold all
    nTone = 1;
    [t1,X1,~] = continuityModel('sustained',nTone,IT(1),0,m,aE,aI,alph,bet,gon,goff,Xinitial);
    [t2,X2,~] = continuityModel('sustained',nTone,IT(2),0,m,aE,aI,alph,bet,gon,goff,Xinitial);
    plot((t1-tStart)/1000,X1,'color',C(1,:),'linewidth',2)
    plot((t2-tStart)/1000,X2,'color',C(2,:),'linewidth',2)
    ylabel('firing rate (x)')
    axis([-.1 1.3 0 1])
    set(gca,'xtick',0:.5:1,'ytick',0:.5:1)
    ti = title('A2');
    set(ti,'position',get(ti,'position')+[-.55 0 0]);
    set(gca,'fontsize',FS)

subplot('position',P(2,:)), hold all
    Inoise = IeqNoise(xeq,aI,alph,INmask);
    idR = find(xeq<=xRnoise(aI,INmask),1,'last'); x1 =xeq(1:idR); I1 = Inoise(1:idR);
    idL = find(xeq>=xLnoise(aI,INmask),1,'first'); x3 =xeq(idL:end); I3 = Inoise(idL:end);
    x2 = xeq(idR:idL); I2 = Inoise(idR:idL);
    hM(1) = plot(I1,x1,'linewidth',3,'color','k');
    hM(2) = plot(I2,x2,':','linewidth',3,'color','k');
    hM(3) = plot(I3,x3,'linewidth',3,'color','k');
    plot(h0(1).XData,h0(1).YData,'linewidth',2,'color',.5*[1 1 1]);
    plot(h0(2).XData,h0(2).YData,':','linewidth',2,'color',.5*[1 1 1]);  
    plot(h0(3).XData,h0(3).YData,'linewidth',2,'color',.5*[1 1 1]);  
    Ieq1 = fzero(@(x) IeqNoise(x,aI,alph,INmask)-IT(1),[1e-8 .1]); plot(IT(1),Ieq1,'o','markeredge','none','markersize',12,'markerfacecolor',C(1,:))
    Ieq2 = fzero(@(x) IeqNoise(x,aI,alph,INmask)-IT(2),[1e-8 .1]); plot(IT(2),Ieq2,'o','markeredge','none','markersize',12,'markerfacecolor',C(2,:))
    axis([-.2 2.2 0 1])
    xlabel('tone strength (I_T)')
    set(gca,'xtick',[0 :1:7],'ytick',0:.5:1,'fontsize',FS)
    ti = title('B1 (masking)'); 
    set(ti,'position',get(ti,'position')+[-.1 0 0]);
    
subplot('position',P(6,:)), hold all
    nTone = 1;
    [t1,X1,~] = continuityModel('sustained',nTone,IT(1),INmask,m,aE,aI,alph,bet,gon,goff,Xinitial);
    [t2,X2,~] = continuityModel('sustained',nTone,IT(2),INmask,m,aE,aI,alph,bet,gon,goff,Xinitial);
    plot((t1-tStart)/1000,X1,'color',C(1,:),'linewidth',2)
    plot((t2-tStart)/1000,X2,'color',C(2,:),'linewidth',2)
    xlabel('time relative to tone onset (sec)')
    axis([-.1 1.3 0 1])
    set(gca,'xtick',0:.5:1,'ytick',0:.5:1)
    ti = title('B2');
    set(ti,'position',get(ti,'position')+[-.55 0 0]);
    set(gca,'fontsize',FS)


subplot('position',P(3,:)), hold all
    Inoise = IeqNoise(xeq,aI,alph,INcont);
    idR = find(xeq<=xRnoise(aI,INcont),1,'last'); x1 =xeq(1:idR); I1 = Inoise(1:idR);
    idL = find(xeq>=xLnoise(aI,INcont),1,'first'); x3 =xeq(idL:end); I3 = Inoise(idL:end);
    x2 = xeq(idR:idL); I2 = Inoise(idR:idL);
    hM(1) = plot(I1,x1,'linewidth',3,'color','k');
    hM(2) = plot(I2,x2,':','linewidth',3,'color','k');
    hM(3) = plot(I3,x3,'linewidth',3,'color','k');
    plot(h0(1).XData,h0(1).YData,'linewidth',2,'color',.5*[1 1 1]);
    plot(h0(2).XData,h0(2).YData,':','linewidth',2,'color',.5*[1 1 1]);  
    plot(h0(3).XData,h0(3).YData,'linewidth',2,'color',.5*[1 1 1]);     
    Ieq1 = fzero(@(x) IeqNoise(x,aI,alph,0)-0,[1e-6 .1]); plot(0,Ieq1,'o','markeredge','none','markersize',12,'markerfacecolor',C(1,:))
    Ieq2 = fzero(@(x) IeqNoise(x,aI,alph,INcont)-0,[.93 .99999]); plot(0,Ieq2,'o','markeredge','none','markersize',12,'markerfacecolor',C(2,:))
    set(gca,'xtick',[0 :1:7],'ytick',0:.5:1,'fontsize',FS)
    ti = title('C1 (continuity)'); 
    set(ti,'position',get(ti,'position')+[+.3 0 0]);
    set(gca,'fontsize',FS)
    axis([-.2 2.2 0 1])
   
subplot('position',P(7,:)), hold all
    nTone = 2;
    [t1,X1,~] = continuityModel('sustained',nTone,IT(2),0,m,aE,aI,alph,bet,gon,goff,Xinitial);
    [t2,X2,~] = continuityModel('sustained',nTone,IT(2),INcont,m,aE,aI,alph,bet,gon,goff,Xinitial);
    plot((t1-tStart)/1000,X1,'color',C(1,:),'linewidth',2)
    i1 = find(t2>=1200,1,'first');
    i2 = find(t2<=1700,1,'last');
    plot((t2(1:i1)-tStart)/1000,X2(1:i1),'color',C(2,:),'linewidth',1)
    plot((t2(i1:i2)-tStart)/1000,X2(i1:i2),'color',C(2,:),'linewidth',2)
    plot((t2(i2:end)-tStart)/1000,X2(i2:end),'color',C(2,:),'linewidth',1)
    axis([-.1 3. 0 1])
    set(gca,'xtick',0:1:5,'ytick',0:.5:1)
    ti = title('C2');
    set(ti,'position',get(ti,'position')+[-1.2 0 0]);
    set(gca,'fontsize',FS)

% threshold curves
aI=1.124; 
alph=.168;

xit = linspace(1,ITmax,20);
for i=1:length(xit)
    maskThresh(i) = fzero( @(IN) xit(i) - IeqNoise(xRnoise(aI,IN), aI, alph, IN) , xit(i)); 
end
contThresh = ones(size(xit))*fzero( @(IN) IeqNoise(xLnoise(aI,IN), aI,alph, IN) , 5);
       
subplot('position',P(4,:)), hold all
    plot(xit,maskThresh,'k:','linewidth',3)
    plot(xit,contThresh,'k','linewidth',3)
    axis([0 ITmax 0 INmax])
    ti = title('D'); ti.Position = ti.Position-[2 0 0];
    set(gca,'fontsize',FS)
    ylabel('I_N at threshold')
    xlabel('I_T')
    text(5.2,7.2,'continuity','fontsize',FS)
    text(5.2,5.4,'masking','fontsize',FS)
    set(gca,'xtick',1:2:5,'ytick',0:5:10)
    
% allowable region
M = fzero( @(IN) ITmax - IeqNoise(xRnoise(aI,IN), aI, alph, IN) , 7); 
C = fzero( @(IN) IeqNoise(xLnoise(aI,IN), aI, alph, IN) , 7);

xai = linspace(.5,3,180);
xalph = linspace(.1,.6,131);linspace(0,1,51);
INC = 5; INM = 5;

for i=1:length(xai)
    for j=1:length(xalph)
        
        maskObj = @(IN) ITmax - IeqNoise(xRnoise(xai(i),IN), xai(i), xalph(j), IN);
        contObj = @(IN) IeqNoise(xLnoise(aI,IN), xai(i), xalph(j), IN);
        
        if (~isreal(maskObj(0)) || ~isreal(maskObj(INmax)))
            MT(i,j) = nan;
        elseif sign(maskObj(0)) == sign(maskObj(INmax))
            MT(i,j) = nan;
        else
            MT(i,j) = fzero( @(IN) maskObj(IN), [realmin INmax]); 
        end

        if (~isreal(contObj(0)) || ~isreal(contObj(INmax)))
            CT(i,j) = nan;
        elseif sign(contObj(0)) == sign(contObj(INmax))
            CT(i,j) = nan;
        else
            CT(i,j) = fzero( @(IN) contObj(IN), [realmin INmax]); 
        end
        if MT(i,j)>0  & MT(i,j)<INmax & CT(i,j)>0 & CT(i,j)<INmax & CT(i,j)>=MT(i,j)
            Z(i,j)= 1;
        else
            Z(i,j) = 0;
        end
    end
end


subplot('position',P(8,:)); hold all;
    W=Z; W(W==0) = nan; MT(isnan(W)) = nan; CT(isnan(W)) = nan; 
    c= contourf(xai,xalph,CT',100,'linestyle','none'); %clabel(c);
    cb = colorbar();
    ylabel(cb,'C(5)')
    cb.Position = [.96 .12 .0115 .3];
    cb.Label.Rotation = 0;
    cb.Label.Position = [1.3 10.9 0];
    c= contour(xai,xalph,MT',[2 3 4 6 8],'k','linewidth',2);% clabel(c);
    text(2.9,.59,'2')
    text(1.91,.395,'3')
    text(1.39,.295,'4')
    text(.93,.21,'6')
    text(.7,.165,'8')
    plot3(aI,alph,100,'ok','markerfacecolor','k','markeredgecolor','k')
    set(gca,'fontsize',FS)
    set(gca,'xtick',1:3,'ytick',.1:.2:.5)
    xlabel('a_I','fontsize',FS)
    ylabel('\alpha','fontsize',FS)
    ti = title('E'); ti.Position = ti.Position-[1 0 0 ];

set(gcf,'units','centimeters','position',[0 0 18 10])


%% MODEL 4: BISTABLE FIGURE

clear all
close all

FS = 10;

figure(4)

dx1 = .23; dx2= .155;
P = [0.07         0.63   dx2    0.3;
     .07+dx1      0.63   dx2    0.3;
     .07+2*dx1    0.63   dx2   0.3;
     .1+3*dx1    0.63   dx2-.045    0.3;
     .07          0.12   dx2    0.3;
     .07+dx1      0.12   dx2    0.3;
     .07+2*dx1    0.12   dx2    0.3;
     .07+3*dx1    0.12   dx2    0.3];
 
 % colors [https://davidmathlogic.com/colorblind/#%23D81B60-%231E88E5-%23FFC107-%23004D40]
 C = [0 90 181; 220 50 32]/255;
 
%%% aE and m calculator
k = 1;
xR = @(aE) 1/2*(1-sqrt(1-4*k./aE)); % right knee at rest
xL = @(aE) 1/2*(1+sqrt(1-4*k./aE)); % lower knee at rest
aEm = @(I) fminsearch(@(P) norm([I(1)-P(2)+log(1/xL(P(1))-1)+P(1)*xL(P(1)) ,  I(2)-P(2)+log(1/xR(P(1))-1)+P(1)*xR(P(1))  ])  , [5 5]);
GON = 5.1856;

IkneeBistable = [-2 2]; % BISTABLE MODEL %
Iknee = IkneeBistable;
a = aEm([Iknee]); aE = a(1); m = a(2);

% equilibrium curve
x=[ linspace(0,.01,1e3), linspace(.01,.99,1e3), linspace(.99,1,1e3)];
Ieq = @(x,m,aE) m - k*log(1./x - 1)-aE*x;

% stimulus parameters
ITmax = 5; INmax = 10;
ITtone = [.8 1.2]; INtone = 0;
ITmask = 3; INmask = [2.8 3.2];
ITcont = 3; INcont = [2.8 3.2];

% model parameters
alph = 0; bet = 2/3; gon = 1; goff = 5.1856; aI = 0;
tStart= 200;
%%%%%% TONE
    nTone = 1;
    Xinitial = [0 , 0];
    [t,X,I] = continuityModel('transient',nTone,ITtone(1),INtone,m,aE,aI,alph,bet,gon,goff,Xinitial);
    t1 = t; X1 = X; I1 = I;
    [t,X,I] = continuityModel('transient',nTone,ITtone(2),INtone,m,aE,aI,alph,bet,gon,goff,Xinitial);
    t2 = t; X2 = X; I2 = I;

subplot('position',P(5,:)); hold all
    plot((t1-tStart)/1000,X1(:,1),'color',C(1,:),'linewidth',2)
    plot((t2-tStart)/1000,X2(:,1),'color',C(2,:),'linewidth',2)
    set(gca,'xtick',0:.5:1)
    ylabel('firing rate (x)')
    axis([-.2 1.2 0 1.])
    ti = title(['A2']);
    ti.Position = [-.100 1.05 0];
    set(gca,'fontsize',FS)

subplot('position',P(1,:)); hold all
    knee1 = fminbnd(@(x) -Ieq(x,m,aE), 0,.5); xEq1 = fzero(@(x) Ieq(x,m,aE),[0.001 knee1]);
    knee2 = fminbnd(@(x) Ieq(x,m,aE),.5,1); xEq3 = fzero(@(x) Ieq(x,m,aE),[knee2 .999]);
    xEq2 = fzero(@(x) Ieq(x,m,aE),[knee1 knee2]);
    
    % stable manifold (linear approx)
    plot(x,aE*(xEq2-x)/GON,'color',.4*[1 1 1],'linewidth',1)
    
    % trajectories in phase plane
    tau = 10; tauEdge = 10;
    i1 = find(t1>=200,1,'first'); i2 = find(t1<1000,1,'last');
    plot(X1(i1:i2,1),(1/GON)*(I1(i1:i2,2)-I1(i1:i2,3)),'color',C(1,:),'linewidth',2)
    i1 = find(t2>=200,1,'first'); i2 = find(t2<1000,1,'last');
    plot(X2(i1:i2,1),(1/GON)*(I2(i1:i2,2)-I2(i1:i2,3)),'color',C(2,:),'linewidth',2)

    % equilibria
    plot(xEq1,0,'o','markeredge','none','markerfacecolor','k','markersize',8)
    plot(xEq2,0,'o','color',.4*[1 1 1],'markersize',8,'linewidth',2)
    plot(xEq3,0,'o','markeredge','none','markerfacecolor','k','markersize',8)

    % vector field
    f     = @(x) 1./(1+exp(-(x-m)/k)); % sigmoidal function
    [xPlane,sScaledPlane] = meshgrid(linspace(0,1,8),linspace(-1.4,1.4,10));
    dx = ( -xPlane + f(aE*xPlane + GON*sScaledPlane) )/tau;
    dsScaled = ( -sScaledPlane)/tauEdge;
    quiver(xPlane,sScaledPlane,dx,dsScaled,.7,'color',.7*[1 1 1])

    % formatting
    ylabel('input variable (s)')
    ti = title('A1 (tone)');
    ti.Position = [.3 1.5 0];
    axis([-.05 1.05 -1.4 1.4])
    set(gca,'fontsize',FS)
   

%%%%%% MASKING
    nTone = 1;
    Xinitial = [0 , 0];
    [t,X,I] = continuityModel('transient',nTone,ITmask,INmask(1),m,aE,aI,alph,bet,gon,goff,Xinitial);
    t1 = t; X1 = X; I1 = I;
    [t,X,I] = continuityModel('transient',nTone,ITmask,INmask(2),m,aE,aI,alph,bet,gon,goff,Xinitial);
    t2 = t; X2 = X; I2 = I;

subplot('position',P(6,:));
    hold all
    plot((t1-tStart)/1000,X1(:,1),'color',C(1,:),'linewidth',2)
    plot((t2-tStart)/1000,X2(:,1),'color',C(2,:),'linewidth',2)
    set(gca,'xtick',0:.5:1)
    xlabel('time relative tone onset (sec)')
    axis([-.2 1.2 0 1.])
    ti = title(['B2']);
    ti.Position = [-.100 1.05 0];
    set(gca,'fontsize',FS)

subplot('position',P(2,:)); 
    hold all
    knee1 = fminbnd(@(x) -Ieq(x,m,aE), 0,.5); xEq1 = fzero(@(x) Ieq(x,m,aE),[0.001 knee1]);
    knee2 = fminbnd(@(x) Ieq(x,m,aE),.5,1); xEq3 = fzero(@(x) Ieq(x,m,aE),[knee2 .999]);
    xEq2 = fzero(@(x) Ieq(x,m,aE),[knee1 knee2]);
    
    % stable manifold (linear approx)
    plot(x,aE*(xEq2-x)/GON,'color',.4*[1 1 1],'linewidth',1)
    
    % trajectories in phase plane
    tau = 10; tauEdge = 10;
    i1 = find(t1>=200,1,'first'); i2 = find(t1<1000,1,'last');
    plot(X1(i1:i2,1),(1/GON)*(I1(i1:i2,2)-I1(i1:i2,3)),'color',C(1,:),'linewidth',2)
    i1 = find(t2>=200,1,'first'); i2 = find(t2<1000,1,'last');
    plot(X2(i1:i2,1),(1/GON)*(I2(i1:i2,2)-I2(i1:i2,3)),'color',C(2,:),'linewidth',2)

    % equilibria
    plot(xEq1,0,'o','markeredge','none','markerfacecolor','k','markersize',8)
    plot(xEq2,0,'o','color',.4*[1 1 1],'markersize',8,'linewidth',2)
    plot(xEq3,0,'o','markeredge','none','markerfacecolor','k','markersize',8)

    % vector field
    f     = @(x) 1./(1+exp(-(x-m)/k)); % sigmoidal function
    [xPlane,sScaledPlane] = meshgrid(linspace(0,1,8),linspace(-1.4,1.4,10));
    dx = ( -xPlane + f(aE*xPlane + GON*sScaledPlane) )/tau;
    dsScaled = ( -sScaledPlane)/tauEdge;
    quiver(xPlane,sScaledPlane,dx,dsScaled,.7,'color',.7*[1 1 1])

    % formatting
    xlabel('firing rate (x)')
    ti = title('B1 (masking)');
    ti.Position = [.4 1.5 0];
    axis([-.05 1.05 -1.4 1.4])
    set(gca,'fontsize',FS)
   


%%%%%% CONTINUITY
    nTone = 2; 
    Xinitial = [0 , 0];
    [t,X,I] = continuityModel('transient',nTone,ITcont,INcont(1),m,aE,aI,alph,bet,gon,goff,Xinitial);
    t1 = t; X1 = X; I1 = I;
    [t,X,I] = continuityModel('transient',nTone,ITcont,INcont(2),m,aE,aI,alph,bet,gon,goff,Xinitial);
    t2 = t; X2 = X; I2 = I;

subplot('position',P(7,:));
    hold all
    plot((t1-tStart)/1000,X1,'color',C(1,:),'linewidth',2)
    i1 = find(t2>=1200,1,'first');
    i2 = find(t2<=1700,1,'last');
    plot((t2(1:i1)-tStart)/1000,X2(1:i1,1),'color',C(2,:),'linewidth',1)
    plot((t2(i1:i2)-tStart)/1000,X2(i1:i2,1),'color',C(2,:),'linewidth',2)
    plot((t2(i2:end)-tStart)/1000,X2(i2:end,1),'color',C(2,:),'linewidth',1)
    axis([-.1 3. 0 1])
    set(gca,'xtick',0:1:5,'ytick',0:.5:1)
    ti = title('C2');
    set(ti,'position',get(ti,'position')+[-1.2 0 0]);
    set(gca,'fontsize',FS)


subplot('position',P(3,:));
    hold all
    knee1 = fminbnd(@(x) -Ieq(x,m,aE), 0,.5); xEq1 = fzero(@(x) Ieq(x,m,aE),[0.001 knee1]);
    knee2 = fminbnd(@(x) Ieq(x,m,aE),.5,1); xEq3 = fzero(@(x) Ieq(x,m,aE),[knee2 .999]);
    xEq2 = fzero(@(x) Ieq(x,m,aE),[knee1 knee2]);
    
    % stable manifold (linear approx)
    plot(x,aE*(xEq2-x)/GON,'color',.4*[1 1 1],'linewidth',1)
    
    % trajectories in phase plane
    tau = 10; tauEdge = 10;
    i1 = find(t1>=1200,1,'first'); i2 = find(t1<1700,1,'last');
    plot(X1(i1:i2,1),(1/GON)*(I1(i1:i2,2)-I1(i1:i2,3)),'color',C(1,:),'linewidth',2)
    i1 = find(t2>=1200,1,'first'); i2 = find(t2<1700,1,'last');
    plot(X2(i1:i2,1),(1/GON)*(I2(i1:i2,2)-I2(i1:i2,3)),'color',C(2,:),'linewidth',2)

    % equilibria
    plot(xEq1,0,'o','markeredge','none','markerfacecolor','k','markersize',8)
    plot(xEq2,0,'o','color',.4*[1 1 1],'markersize',8,'linewidth',2)
    plot(xEq3,0,'o','markeredge','none','markerfacecolor','k','markersize',8)

    % vector field
    f     = @(x) 1./(1+exp(-(x-m)/k)); % sigmoidal function
    [xPlane,sScaledPlane] = meshgrid(linspace(0,1,8),linspace(-1.4,1.4,10));
    dx = ( -xPlane + f(aE*xPlane + GON*sScaledPlane) )/tau;
    dsScaled = ( -sScaledPlane)/tauEdge;
    quiver(xPlane,sScaledPlane,dx,dsScaled,.7,'color',.7*[1 1 1])

    % formatting
    ti = title('C1 (continuity)');
    ti.Position = [.5 1.5 0];
    axis([-.05 1.05 -1.4 1.4])
    set(gca,'fontsize',FS)
   
    
% THRESHOLDS

xIT = linspace(1,5,3);
subplot('position',P(4,:)); 
    hold all
    plot(xIT, (1/bet)*(xIT-1),'k','linewidth',3)
    axis([0 ITmax 0 INmax])
    ti = title('D'); ti.Position = ti.Position-[2 0 0];
    set(gca,'fontsize',FS)
    ylabel('I_N at threshold')
    xlabel('I_T')
    set(gca,'xtick',1:2:5,'ytick',0:5:10)
    text(5.2,7.3,'continuity')
    text(5.2,5.3,'= masking')
    axis([0 ITmax 0 INmax])
    set(gca,'fontsize',FS)

    
    % LEGEND
    annotation('textbox',[.75 .35 .1 .1],'String','legend','linestyle','none','fontsize',FS,'fontweight','bold')

    annotation('textbox',[.75 .285 .1 .1],'String','tone only (A)','linestyle','none','fontsize',FS)
    annotation('line',[.76 .8],[.30 .30],'linewidth',2,'color',C(1,:))
    annotation('textbox',[.8 .235 .1 .1],'String',['I_T = ',num2str(ITtone(1)), ', I_N = 0'],'linestyle','none','fontsize',FS)
    annotation('line',[.76 .8],[.255 .255],'linewidth',2,'color',C(2,:))
    annotation('textbox',[.8 .185 .1 .1],'String',['I_T = ',num2str(ITtone(2)), ', I_N = 0'],'linestyle','none','fontsize',FS)

    annotation('textbox',[.75 .11 .1 .1],'String','tone with noise (B & C)','linestyle','none','fontsize',FS)
    annotation('line',[.76 .8],[.13 .13],'linewidth',2,'color',C(1,:))
    annotation('textbox',[.8 .07 .1 .1],'String',['I_T = ',num2str(ITmask), ', I_N = ',num2str(INmask(1))],'linestyle','none','fontsize',FS)
    annotation('line',[.76 .8],[.08 .08],'linewidth',2,'color',C(2,:))
    annotation('textbox',[.8 .01 .1 .1],'String',['I_T = ',num2str(ITmask), ', I_N = ',num2str(INmask(2))],'linestyle','none','fontsize',FS)


set(gcf,'units','centimeters','position',[0 0 18 10])
set(gcf, 'PaperPositionMode','auto')




%% FIGURE 5: combined inputs demo



clear all
close all; 

figure(5)

P = [0.08    0.21    0.36    0.7;
     0.6    0.21    0.36    0.7];

%%% aE and m calculator
k = 1;
xR = @(aE) 1/2*(1-sqrt(1-4*k./aE)); % right knee at rest
xL = @(aE) 1/2*(1+sqrt(1-4*k./aE)); % lower knee at rest
aEm = @(I) fminsearch(@(P) norm([I(1)-P(2)+log(1/xL(P(1))-1)+P(1)*xL(P(1)) ,  I(2)-P(2)+log(1/xR(P(1))-1)+P(1)*xR(P(1))  ])  , [5 5]);

IkneeCombined = [.2 6]; % COMBINED MODEL %
Iknee = IkneeCombined;
a = aEm([Iknee]); aE = a(1); m = a(2);

% equilibrium curve
x=[ linspace(0,.01,1e3), linspace(.01,.99,1e3), linspace(.99,1,1e3)];
Ieq = @(x,m,aE) m - k*log(1./x - 1)-aE*x;

tStart = 200;
C = lines(3);
subplot('position',P(1,:)); hold all

    IT = 1.5;
    I = Ieq(x,m,aE);
    idR = find(x<=xR(aE),1,'last'); x1 =x(1:idR); I1 = I(1:idR);
    idL = find(x>=xL(aE),1,'first'); x3 =x(idL:end); I3 = I(idL:end);
    x2 = x(idR:idL); I2 = I(idR:idL);
    h(1) = plot(I1,x1,'linewidth',3,'color','k');
    h(2) = plot(I2,x2,':','linewidth',3,'color','k');
    h(3) = plot(I3,x3,'linewidth',3,'color','k');
    I0 = fzero(@(x) Ieq(x,m,aE),0.0005); plot(0,I0,'o','markeredge','none','markersize',12,'markerfacecolor',C(2,:))%'markerfacecolor','k')
    I2 = fzero(@(x) Ieq(x,m,aE)-IT,[1e-6,.1]); plot(IT,I2,'o','markeredge','none','markersize',12,'markerfacecolor',C(1,:))
    I3 = fzero(@(x) Ieq(x,m,aE)-IT,[.9 .9999]); plot(IT,I3,'o','markeredge','none','markersize',12,'markerfacecolor',C(3,:))
    xlabel('tone level (I_T)')
    ylabel('firing rate (x)')
    axis([-.3 5 0 1.1])
    set(gca,'xtick',[0 :1:5],'ytick',0:.2:1)
    ti = title(['A']);
    ti.Position = [-1.1 1.1 0];

    CC = [0 0 0 ; C ; 0 0 0];
    subplot('position',P(2,:)); hold all
    alph = 0; bet = 0; aI = 0; IN = 0; nTone = 1; Xinit = [0 0];

    [t,X,I] = TwoPopFunction('transient0',nTone,IT,IN,m,aE,aI,alph,bet,1,[],Xinit);
    plot((t-tStart)/1000,X(:,1),'linewidth',2,'color',C(2,:))

    [t,X,I] = TwoPopFunction('combined',nTone,IT,IN,m,aE,aI,alph,bet,1,[],Xinit);
    plot((t-tStart)/1000,X(:,1),'linewidth',2,'color',C(3,:))

    [t,X,I] = TwoPopFunction('sustained',nTone,IT,IN,m,aE,aI,alph,bet,0,0,Xinit);
    plot((t-tStart)/1000,X(:,1),'linewidth',2,'color',C(1,:))

    xlabel('time relative to tone onset (s)')
    ylabel('firing rate (x)')
    axis([-.2 1.300 0 1.1])
    ti = title(['B']);
    ti.Position = [-.230 1.1 0];

set(gcf,'units','centimeters','position',[0 0 18 6])
set(gcf, 'PaperPositionMode','auto') 






%% FIGURE 6: combined figure

clear all
close all

FS = 10;

figure(6)

dx1 = .23; dx2= .155;
P = [0.08    0.63   dx2    0.3;
     .08+dx1    0.63   dx2    0.3;
     .08+2*dx1    0.63   dx2   0.3;
     .1+3*dx1    0.63   dx2-.04    0.3;
     .08    0.12   dx2    0.3;
     .08+dx1    0.12   dx2    0.3;
     .08+2*dx1    0.12   dx2    0.3;
     .1+3*dx1    0.12   dx2    0.3];
 
 % colors [https://davidmathlogic.com/colorblind/#%23D81B60-%231E88E5-%23FFC107-%23004D40]
 C = [0 90 181; 220 50 32]/255;

%%% aE and m calculator
k = 1;
xR = @(aE) 1/2*(1-sqrt(1-4*k./aE)); % right knee at rest
xL = @(aE) 1/2*(1+sqrt(1-4*k./aE)); % lower knee at rest
aEm = @(I) fminsearch(@(P) norm([I(1)-P(2)+log(1/xL(P(1))-1)+P(1)*xL(P(1)) ,  I(2)-P(2)+log(1/xR(P(1))-1)+P(1)*xR(P(1))  ])  , [5 5]);


IkneeCombined = [.2 6]; % COMBINED MODEL %
Iknee = IkneeCombined;
a = aEm([Iknee]); aE = a(1); m = a(2);

% equilibrium curve
x=[ linspace(0,.00001,1e3), linspace(.00001,.99999,1e3), linspace(.99999,1,1e3)];
f   = @(x) 1./(1+exp(-(x-m)/k)); % sigmoidal function
Ieq = @(x,m,aE) m - k*log(1./x - 1)-aE*x;
IeqNoise = @(x,aI,alph,IN) m - k*log(1./x - 1)-aE*x +  aI*IN*(1-x) -alph*IN;


% gOn calcualator
I0 =1; % for combined. threshold at sustained input = 1
[xMin,~] = fminbnd(@(x) (-x + f(aE*x+I0)), 0,1); 
[xMax,~] = fminbnd(@(x) -(-x + f(aE*x+I0)), 0,1); 
xR = fzero( @(x) -x + f(aE*x+I0), [1e-4,xMin-1e-4]); % rest equilibrium
xS = fzero( @(x) -x + f(aE*x+I0), [xMin+1e-4,xMax-1e-4]); % saddle equilibrium
gon = aE*(xS-xR);

ITmax = 5; INmax = 10;
maskTarget = 6;
contTarget = 7;

tStart = 200;

% minimizing forward maskingregion
alph = .5; aI = 7;
bet = fzero(@(bet) combinedMaskThresholdFunction(ITmax,maskTarget,m,aE,aI,alph,bet,gon,0), [0 1]);
goff =fzero( @(goff) combinedContThresholdFunction(ITmax,contTarget,m,aE,aI,alph,bet,gon,goff), [0 20]);


%%%%%% TONE
    IN = 0; IT = [1.2]; nTone = 1;
    Xinitial = [0 , 0];
    [t,X,I] = continuityModel('combined',nTone,IT,0,m,aE,aI,alph,bet,gon,goff,Xinitial);
    t1 = t; X1 = X; I1 = I;

subplot('position',P(5,:)); hold all
    plot((t1-tStart)/1000,X1(:,1),'color',C(1,:),'linewidth',2)
    set(gca,'xtick',0:.5:1)
    ylabel('firing rate (x)')
    axis([-.2 1.2 0 1.])
    text(.25,.29,['I_T = ',num2str(IT)])
    text(.25,.12,['I_N = ',num2str(0)])
    ti = title(['A2']);
    ti.Position = [-.1 1. 0];
    set(gca,'fontsize',FS)

subplot('position',P(1,:)); hold all
    
    knee1 = fminbnd(@(x) -Ieq(x,m,aE)+IT, 0,.5); xEq1 = fzero(@(x) Ieq(x,m,aE)-IT,[0.00001 knee1]);
    knee2 = fminbnd(@(x) Ieq(x,m,aE)-IT,.5,1); xEq3 = fzero(@(x) Ieq(x,m,aE)-IT,[knee2 .99999]);
    xEq2 = fzero(@(x) Ieq(x,m,aE)-IT,[knee1 knee2]);
    
    % stable manifold (linear approx)
    xEqInactive0 = fzero(@(x) Ieq(x,m,aE)-0,[0.000001 .999]);
    xEqSaddle1 = fzero(@(x) Ieq(x,m,aE)-1,[0.1 .9]);
    plot(x,(xEq2-x)/(xEqSaddle1-xEqInactive0),'color',.4*[1 1 1],'linewidth',1)
    
    % trajectories in phase plane
    tau = 10; tauEdge = 10;
    i1 = find(t1>=200,1,'first'); i2 = find(t1<1000,1,'last');
    plot(X1(i1:i2,1),(1/gon)*(I1(i1:i2,2)-I1(i1:i2,3)),'color',C(1,:),'linewidth',2)

    % equilibria
    plot(xEq1,0,'o','markeredge','none','markerfacecolor','k','markersize',8)
    plot(xEq2,0,'o','color',.4*[1 1 1],'markersize',8,'linewidth',2)
    plot(xEq3,0,'o','markeredge','none','markerfacecolor','k','markersize',8)

    % vector field
    f     = @(x) 1./(1+exp(-(x-m)/k)); % sigmoidal function
    [xPlane,sScaledPlane] = meshgrid(linspace(0,1,8),linspace(-1.4,1.4,10));
    dx = ( -xPlane + f(aE*xPlane + gon*sScaledPlane + IT) )/tau;
    dsScaled = ( -sScaledPlane)/tauEdge;
    quiver(xPlane,sScaledPlane,dx,dsScaled,.7,'color',.7*[1 1 1])

    % formatting
    ylabel('input variable (s)')
    ti = title('A1 (tone)');
    ti.Position = [.3 1.5 0];
    axis([-.05 1.05 -1.4 1.4])
    set(gca,'fontsize',FS)
   
    
% MASKING   
    IN = [1.5]; IT = 2;
    nTone = 1;
    Xinitial = [0 , 0];
    [t,X,I] = continuityModel('combined',nTone,IT,IN,m,aE,aI,alph,bet,gon,goff,Xinitial);
    t1 = t; X1 = X; I1 = I;

subplot('position',P(6,:)); hold all
    plot((t1-tStart)/1000,X1(:,1),'color',C(1,:),'linewidth',2)
    xlabel('time relative tone onset (sec)')
    axis([-.2 1.2 0 1.])
    text(.35,.29,['I_T = ',num2str(IT)])
    text(.35,.12,['I_N = ',num2str(IN)])
    ti = title(['B2']);    
    ti.Position = [-.1 1. 0];
    set(gca,'fontsize',FS)
    
subplot('position',P(2,:)); hold all

    % equilibria
    cut1 = find((IeqNoise(x,aI,alph,IN)-IT)>0,1,'first');
    cut2 = find((IeqNoise(x,aI,alph,IN)-IT)<0,1,'last');
    xEqSaddle  = fzero( @(x) IeqNoise(x,aI,alph,IN)-IT, [x(cut1) x(cut2)]);
    xEqI = fzero( @(x) IeqNoise(x,aI,alph,IN)-IT, [1e-12 x(cut1)*(1-eps) ]);
    xEqA = fzero( @(x) IeqNoise(x,aI,alph,IN)-IT, [x(cut2)*(1+eps) .999999]);
    plot(xEqI,0,'o','markeredge','none','markerfacecolor','k','markersize',8)
    plot(xEqSaddle,0,'o','color',.4*[1 1 1],'markersize',8,'linewidth',2)
    plot(xEqA,0,'o','markeredge','none','markerfacecolor','k','markersize',8)
    
    % stable manifold (linear approx)
    cut1 = find((IeqNoise(x,aI,alph,0)-0)>0,1,'first');
    cut2 = find((IeqNoise(x,aI,alph,0)-0)<0,1,'last');
    xEqInactive0  = fzero( @(x) IeqNoise(x,aI,alph,0)-0, [x(cut1) x(cut2)]);

    cut1 = find((IeqNoise(x,aI,alph,0)-1)>0,1,'first');
    cut2 = find((IeqNoise(x,aI,alph,0)-1)<0,1,'last');
    xEqSaddle1  = fzero( @(x) IeqNoise(x,aI,alph,0)-1, [x(cut1) x(cut2)]);
    
    cut1 = find((IeqNoise(x,aI,alph,IN)-IT)>0,1,'first');
    cut2 = find((IeqNoise(x,aI,alph,IN)-IT)<0,1,'last');
    xEqSaddle  = fzero( @(x) IeqNoise(x,aI,alph,IN)-IT, [x(cut1) x(cut2)]);
    plot(x,(1+aI*IN/aE)*(xEqSaddle-x)/(xEqSaddle1-xEqInactive0),'color',.4*[1 1 1],'linewidth',1)

    % trajectories in phase plane
    tau = 10; tauEdge = 10;
    i1 = find(t1>=200,1,'first'); i2 = find(t1<1000,1,'last');
    plot(X1(i1:i2,1),(1/gon)*(I1(i1:i2,2)-I1(i1:i2,3)),'color',C(1,:),'linewidth',2)

    % vector field
    f     = @(x) 1./(1+exp(-(x-m)/k)); % sigmoidal function
    [xPlane,sScaledPlane] = meshgrid(linspace(0,1,8),linspace(-1.4,2.2,10));
    dx = ( -xPlane + f(aE*xPlane -aI*IN*(1-xPlane) + gon*sScaledPlane + IT + alph*IN) )/tau;
    dsScaled = ( -sScaledPlane)/tauEdge;
    quiver(xPlane,sScaledPlane,dx,dsScaled,.7,'color',.7*[1 1 1])

    % formatting
    xlabel('firing rate (x)')
    ti = title('B1 (masking)');
    ti.Position = [.4 2.3 0];
    axis([-.05 1.05 -1.4 2.2])
    set(gca,'fontsize',FS)
        
%%%%%%% CONTINUITY
    IN = [3]; IT = 2;
    nTone = 2;
    Xinitial = [0 , 0];
    [t,X,I] = continuityModel('combined',nTone,IT,IN,m,aE,aI,alph,bet,gon,goff,Xinitial);
    t1 = t; X1 = X; I1 = I;

subplot('position',P(7,:)); hold all
    plot((t1-tStart)/1000,X1(:,1),'color',C(1,:),'linewidth',2)
    text(.8,.29,['I_T = ',num2str(IT)])
    text(.8,.12,['I_N = ',num2str(IN)])
    axis([-.1 3. 0 1])
    set(gca,'xtick',0:1:5,'ytick',0:.5:1)
    ti = title('C2');
    set(ti,'position',get(ti,'position')+[-1.2 0 0]);
    set(gca,'fontsize',FS)

    
subplot('position',P(3,:)); hold all
    % equilibria
    cut1 = find((IeqNoise(x,aI,alph,IN)-0)>0,1,'first');
    cut2 = find((IeqNoise(x,aI,alph,IN)-0)<0,1,'last');
    xEqSaddle  = fzero( @(x) IeqNoise(x,aI,alph,IN)-0, [x(cut1) x(cut2)]);
    xEqI = fzero( @(x) IeqNoise(x,aI,alph,IN)-0, [realmin x(cut1)*(1-eps) ]);
    xEqA = fzero( @(x) IeqNoise(x,aI,alph,IN)-0, [x(cut2)*(1+eps) 1-eps]);
    plot(xEqI,0,'o','markeredge','none','markerfacecolor','k','markersize',8)
    plot(xEqSaddle,0,'o','color',.4*[1 1 1],'markersize',8,'linewidth',2)
    plot(xEqA,0,'o','markeredge','none','markerfacecolor','k','markersize',8)
    
    % stable manifold (linear approx)    
    cut1 = find((IeqNoise(x,aI,alph,IN)-0)>0,1,'first');
    cut2 = find((IeqNoise(x,aI,alph,IN)-0)<0,1,'last');
    xEqSaddle  = fzero( @(x) IeqNoise(x,aI,alph,IN)-0, [x(cut1) x(cut2)]);
    plot(x,(aE+aI*IN)*(xEqSaddle-x)/goff,'color',.4*[1 1 1],'linewidth',1)

    % trajectories in phase plane
    tau = 10; tauEdge = 10;
    i1 = find(t1>=1200,1,'first'); i2 = find(t1<1700,1,'last');
    plot(X1(i1:i2,1),(1/goff)*(I1(i1:i2,2)-I1(i1:i2,3)),'color',C(1,:),'linewidth',2)

    % vector field
    f     = @(x) 1./(1+exp(-(x-m)/k)); % sigmoidal function
    [xPlane,sScaledPlane] = meshgrid(linspace(0,1,8),linspace(-2.4,1.4,10));
    dx = ( -xPlane + f(aE*xPlane -aI*IN*(1-xPlane) + goff*sScaledPlane + alph*IN) )/tau;
    dsScaled = ( -sScaledPlane)/tauEdge;
    quiver(xPlane,sScaledPlane,dx,dsScaled,.7,'color',.7*[1 1 1])
 
    % formatting
    ti = title('C1 (continuity)');
    ti.Position = [.5 1.5 0];
    axis([-.05 1.05 -2.4 1.4])
    set(gca,'fontsize',FS)

    
% THRESHOLDS

xIT = linspace(1,5);
subplot('position',P(4,:)), 
hold all
for i=1:length(xIT)
    cIN(i) = fzero(@(IN) combinedContThresholdFunction(xIT(i),IN,m,aE,aI,alph,bet,gon,goff), [0 INmax]);
    mIN(i) = fzero(@(IN) combinedMaskThresholdFunction(xIT(i),IN,m,aE,aI,alph,bet,gon,goff), [0 INmax]);
    fIN(i) = fzero(@(IN) combinedForwardMaskThresholdFunction(xIT(i),IN,m,aE,aI,alph,bet,gon,goff), [0 200]);
end
    plot(xIT,mIN,'k:','linewidth',2)
    plot(xIT,cIN,'k','linewidth',2)
    axis([0 ITmax 0 INmax])
    ti = title('D'); ti.Position = ti.Position-[2 0 0];
    set(gca,'fontsize',FS)
    ylabel('I_N at threshold')
    xlabel('I_T')
    text(5.2,7.4,'continuity','fontsize',FS)
    text(5.2,5.8,'masking','fontsize',FS)
    set(gca,'xtick',1:2:5,'ytick',0:5:10)
    
    


subplot('position',P(8,:)), 
hold all
for i=1:length(xIT)
    cINlow(i) = fzero(@(IN) combinedContThresholdFunction(xIT(i),IN,m,aE,aI,alph,bet,gon,goff*.75), [0 10]);
    cINhigh(i) = fzero(@(IN) combinedContThresholdFunction(xIT(i),IN,m,aE,aI,alph,bet,gon,goff*1.25), [0 20]);
end
plot(xIT, cIN,'color','k','linewidth',2)
plot(xIT, cINlow,'color',C(1,:),'linewidth',2)
plot(xIT, cINhigh,'color',C(2,:),'linewidth',2)
    axis([0 ITmax 0 INmax])
    ti = title('E'); ti.Position = ti.Position-[2 0 0];
    set(gca,'fontsize',FS)
    ylabel('C(I_N)')
    xlabel('I_T')
    text(5.2,10.7,'\gamma_{off}','fontsize',FS)
    text(5.2,9,'1.10','fontsize',FS)
    text(5.2,7,'0.88','fontsize',FS)
    text(5.2,5,'0.66','fontsize',FS)
    set(gca,'xtick',1:2:5,'ytick',0:5:10)
    set(gca,'fontsize',FS)

set(gcf,'units','centimeters','position',[0 0 18 10])
set(gcf, 'PaperPositionMode','auto') 




%% USEFUL FUNCTIONS


function out = combinedForwardMaskThresholdFunction(IT,IN,m,aE,aI,alph,bet,gon,goff)
    
    ITmax = 5;
    INmax = 10 ;

    k=1;
    tau = 10; tauEdge = 10;

    % sigmoidal function
    f = @(x) 1./(1+exp(-(x-m)/k)); 

    % x-I equilibrium curve, with noise
    fI = @(x,IT,IN) f(aE*x + IT + alph*IN - aI*IN*(1-x));

    % find inactive state: no input
    xRough = 1e-6:.01:1-1e-6;
    [~,a] = min(-xRough + fI(xRough,0,0)); xMin = xRough(a);
    xI0 = fzero( @(x) -x +  fI(x,0,0), [0,xMin-1e-4]); % rest equilibrium

    % find saddle: tone=1
    [~,a] = min(-xRough + fI(xRough,1,0)); xMin = xRough(a);
    [~,b] = max(-xRough + fI(xRough,1,0)); xMax = xRough(b);
    xS1 = fzero( @(x) -x +  fI(x,1,0), [xMin+1e-4,xMax-1e-4]); % saddle equilibrium
    
    % find saddle: tone only
    [~,a] = min(-xRough + fI(xRough,IT,0)); xMin = xRough(a);
    [~,b] = max(-xRough + fI(xRough,IT,0)); xMax = xRough(b);
    xS = fzero( @(x) -x +  fI(x,IT,0), [xMin+1e-4,xMax-1e-4]); % saddle equilibrium

    % for forward masking threshold. activation from noise only inactive
    [~,a] = min(-xRough + fI(xRough,0,IN)); xMin = xRough(a);
    xMask = fzero( @(x) -x +  fI(x,0,IN), [0,xMin-1e-4]); % rest equilibrium
    
    % separatrix set up for root finding
    out = max(IT-bet*IN,0) - ( (1)*(xS - xMask) / (xS1-xI0) ) ;
    
end



function out = combinedMaskThresholdFunction(IT,IN,m,aE,aI,alph,bet,gon,goff)
    
    ITmax = 5;
    INmax = 10 ;

    k=1;
    tau = 10; tauEdge = 10;

    % sigmoidal function
    f = @(x) 1./(1+exp(-(x-m)/k)); 

    % x-I equilibrium curve, with noise
    fI = @(x,IT,IN) f(aE*x + IT + alph*IN - aI*IN*(1-x));

    % find inactive state: no input
    xRough = 1e-6:.01:1-1e-6;
    [~,a] = min(-xRough + fI(xRough,0,0)); xMin = xRough(a);
    xI0 = fzero( @(x) -x +  fI(x,0,0), [0,xMin-1e-4]); % rest equilibrium

    % find saddle: tone=1
    [~,a] = min(-xRough + fI(xRough,1,0)); xMin = xRough(a);
    [~,b] = max(-xRough + fI(xRough,1,0)); xMax = xRough(b);
    xS1 = fzero( @(x) -x +  fI(x,1,0), [xMin+1e-4,xMax-1e-4]); % saddle equilibrium
    
    % find saddle: tone and noise
    [~,a] = min(-xRough + fI(xRough,IT,IN)); xMin = xRough(a);
    [~,b] = max(-xRough + fI(xRough,IT,IN)); xMax = xRough(b);
    xS = fzero( @(x) -x +  fI(x,IT,IN), [xMin+1e-4,xMax-1e-4]); % saddle equilibrium

    % for masking threshold. activation from rest
    xMask = xI0;
    
    % separatrix set up for root finding
    out = max(IT-bet*IN,0) - ( (1 + aI*IN/aE)*(xS - xMask) / (xS1-xI0) ) ;
    
end




function out = combinedContThresholdFunction(IT,IN,m,aE,aI,alph,bet,gon,goff)

    ITmax = 5;
    INmax = 10 ;

    k=1;
    tau = 10; tauEdge = 10;

    % sigmoidal function
    f = @(x) 1./(1+exp(-(x-m)/k)); 

    % x-I equilibrium curve, with noise
    fI = @(x,IT,IN) f(aE*x + IT + alph*IN - aI*IN*(1-x));

    % find inactive state: no input
    xRough = linspace(realmin,1-eps,100);
    
    % find saddle: continuity: noise only
    [~,a] = min(-xRough + fI(xRough,0,IN)); xMin = xRough(a);
    [~,b] = max(-xRough + fI(xRough,0,IN)); xMax = xRough(b);
    if length(find(diff(sign(-xRough + fI(xRough,0,IN)))))<2;
        out = 999; % no sadddle
    else
        xS = fzero( @(x) -x +  fI(x,0,IN), [xMin,xMax]); % saddle equilibrium

        % for continuity: deactivation from active, tone only
        [~,a] = max(-xRough + fI(xRough,IT,0)); xMax = xRough(a);
        xCont = fzero( @(x) -x +  fI(x,IT,0), [xMax+1e-4,1]); % active equilibrium

        % separatrix set up for root finding
        out = goff*max(IT-bet*IN,0) + (aE + aI*IN)*(xS - xCont)   ;
    end

end

