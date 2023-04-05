%CellularGrowthFluctNutrientEnvF:
ColorSeq=[0 .45 .74 ; .85 .33 .1; .93 .69 .13; .49 .18 .56; .47 .67 .19];

%% I. Memory capacity effect on transition time and estimation accuracy
rL=.05;
rH=1;
p=.5;
betaCrit=1+rH/rL
beta=1+rL
pI=(beta-1)*rL/((beta-1)*rL+rH)
N=1000;
method=NaN;
Hx=1;
NTraj=10;
lMax=inf;
lMin=1;

NIterate=5*10^3;
LSet=[1 2 5 25 50]

StoreRAvgFixedLConstantToFluctuating=cell(1,length(LSet));
Hits=nan(1,length(LSet));
t=1;
for l=LSet
    for z=1:NTraj
        [Env R]=StochasticSensing1(rL, rH, beta, p, l, N, method, Hx);
    end
    %Characterize average response
    RTotals  =nan(NIterate,N+l);
    HitTotals=nan(NIterate,N);
    parfor z=1:NIterate
        [Env R Hit]=StochasticSensing1(rL, rH, beta, p, l, N, method, Hx);
        RTotals(z,:)=R;
        HitTotals(z,:)=Hit;
    end
    RAvg=mean(RTotals,1);
    StoreRAvgFixedLConstantToFluctuating{t}=RAvg;

    kCrit(t)=pI/p*l;

    pMiss=p*(1-p)/((p*(1-p)+l*(p-pI)));
    RF=rL*(1-p)+rH*p;
    RL=beta*rL*(1-p);
    Rlower(t)=( rH*p-(beta-1)*rL*(1-p) )* pMiss;
    Rlower2(t)= RF - (RF+RL)*pMiss;
    
    CDFBinom=binocdf(pI*l,l,p);
    Rlower3(t)= RF*(1-CDFBinom)+RL*CDFBinom;
    
    t=t+1;
end

figure; hold on; box on;
for z=1:length(LSet)
    RAvg=StoreRAvgFixedLConstantToFluctuating{z};
    plot([-LSet(z)+1:1:N],RAvg,'LineWidth',2);
end
plot([1:1:N],(rL+rH)/2*ones(N,1),'k-','LineWidth',1.5); hold on;
plot(0,0,'k--','LineWidth',1.5);
plot(0,0,'k*','Linewidth',1.5);
plot([1:1:N],rL*ones(N,1),'k-','LineWidth',1.5);

for z=1:length(LSet)
    plot(max(1,kCrit(z)),rL,'*','Color',ColorSeq(z,:),'LineWidth',1.5);
end
for z=1:length(LSet)
    PredDiff3=Rlower3(z);
    plot([1 N],PredDiff3*[1 1],'--','Color',.5*ColorSeq(z,:),'LineWidth',1.5);
end
xlim([1 N]); set(gca,'XScale','log');
hl=legend(sprintf('$\\ell=$%.f',LSet(1)),sprintf('$\\ell=$%.f',LSet(2)),...
    sprintf('$\\ell=$%.f',LSet(3)), sprintf('$\\ell=$%.f',LSet(4)),...
    sprintf('$\\ell=$%.f',LSet(5)),'limits','theory prediciton','inflection','location','Southeast');
set(hl,'interpreter','latex');
xlabel('Time','interpreter','latex');
ylabel('Mean growth rate','interpreter','latex');
title('Growth: Constant-to-fluctuating environments, $\beta=1+r_L$','interpreter','latex');

%2. beta=beta_crit/5
rL=.05;
rH=1;
p=.5;
betaCrit=1+rH/rL
beta=betaCrit/5
pI=(beta-1)*rL/((beta-1)*rL+rH);
N=10^3;
l=5;
method=NaN;
Hx=1;
NTraj=10;
lMax=inf;
lMin=1;

NIterate=5*10^3;

LSet=[1 2 5 25 50]
StoreRAvgFixedLConstantToFluctuating=cell(1,length(LSet));
Hits=nan(1,length(LSet));
t=1;
tic;
for l=LSet
    toc; l
    for z=1:NTraj
        [Env R]=StochasticSensing1(rL, rH, beta, p, l, N, method, Hx);
    end
   
    RTotals  =nan(NIterate,N+l);
    HitTotals=nan(NIterate,N);
    parfor z=1:NIterate
        [Env R Hit]=StochasticSensing1(rL, rH, beta, p, l, N, method, Hx);
        RTotals(z,:)=R;
        HitTotals(z,:)=Hit;
    end
    RAvg=mean(RTotals,1);
    StoreRAvgFixedLConstantToFluctuating{t}=RAvg;

    kCrit(t)=pI/p*l;

    pMiss=p*(1-p)/((p*(1-p)+l*(p-pI)));
    RF=rL*(1-p)+rH*p;
    RL=beta*rL*(1-p);
    Rlower(t)=( rH*p-(beta-1)*rL*(1-p) )* pMiss;
    Rlower2(t)=RF - (RF+RL)*pMiss;
    
    CDFBinom=binocdf(pI*l,l,p);
    Rlower3(t)= RF*(1-CDFBinom)+RL*CDFBinom;
    
    t=t+1;
end


figure; hold on; box on;
for z=1:length(LSet)
    RAvg=StoreRAvgFixedLConstantToFluctuating{z};
    plot([-LSet(z)+1:1:N],RAvg,'LineWidth',2);
end
plot([1:1:N],(rL+rH)/2*ones(N,1),'k-','LineWidth',1.5); hold on;
plot(0,0,'k--','LineWidth',1.5);
plot(0,0,'k*','Linewidth',1.5);

plot([1:1:N],rL*ones(N,1),'k-','LineWidth',1.5);

for z=1:length(LSet)
    plot(max(1,kCrit(z)),rL,'*','Color',ColorSeq(z,:),'LineWidth',1.5);
end
for z=1:length(LSet)
    PredDiff3=Rlower3(z);
    plot([1 N],PredDiff3*[1 1],'--','Color',.5*ColorSeq(z,:),'LineWidth',1.5);
end
xlim([1 N]); ylim([0 .6]);
set(gca,'XScale','log');
hl=legend(sprintf('$\\ell=$%.f',LSet(1)),sprintf('$\\ell=$%.f',LSet(2)),...
    sprintf('$\\ell=$%.f',LSet(3)), sprintf('$\\ell=$%.f',LSet(4)),...
    sprintf('$\\ell=$%.f',LSet(5)),'limits','theory prediction','inflection','location','Southeast');
set(hl,'interpreter','latex');
xlabel('Time','interpreter','latex');
ylabel('Mean growth rate','interpreter','latex');
title('Growth: Constant-to-fluctuating environments, $\beta=\beta_{Crit}/5$','interpreter','latex');


%3. beta=betaCrit/2
rL=.05;
rH=1;
p=.5;
betaCrit=1+rH/rL
beta=betaCrit/2
pI=(beta-1)*rL/((beta-1)*rL+rH);
N=1000;
l=5;
method=NaN;
Hx=1;
NTraj=10;
lMax=inf;
lMin=1;

NIterate=5*10^3;
LSet=[1 2 5 25 50]
StoreRAvgFixedLConstantToFluctuating=cell(1,length(LSet));
Hits=nan(1,length(LSet));
t=1;
for l=LSet
    for z=1:NTraj
        [Env R]=StochasticSensing1(rL, rH, beta, p, l, N, method, Hx);
    end

    RTotals  =nan(NIterate,N+l);
    HitTotals=nan(NIterate,N);
    parfor z=1:NIterate
        [Env R Hit]=StochasticSensing1(rL, rH, beta, p, l, N, method, Hx);
        RTotals(z,:)=R;
        HitTotals(z,:)=Hit;
    end
    RAvg=mean(RTotals,1);
    StoreRAvgFixedLConstantToFluctuating{t}=RAvg;

    kCrit(t)=pI/p*l;
    
    pMiss=p*(1-p)/((p*(1-p)+l*(p-pI)));
    RF=rL*(1-p)+rH*p;
    RL=beta*rL*(1-p);
    Rlower(t)=( rH*p-(beta-1)*rL*(1-p) )* pMiss;
    Rlower2(t)=RF - (RF+RL)*pMiss;
    
    CDFBinom=binocdf(pI*l,l,p);
    Rlower3(t)= RF*(1-CDFBinom)+RL*CDFBinom;
    
    t=t+1;
end


figure; hold on; box on;
for z=1:length(LSet)
    RAvg=StoreRAvgFixedLConstantToFluctuating{z};
    plot([-LSet(z)+1:1:N],RAvg,'LineWidth',2);
end
plot([1:1:N],(rL+rH)/2*ones(N,1),'k-','LineWidth',1.5); hold on;
plot(0,0,'k--','LineWidth',1.5);
plot(0,0,'k*','Linewidth',1.5);

plot([1:1:N],rL*ones(N,1),'k-','LineWidth',1.5);

for z=1:length(LSet)
    plot(max(1,kCrit(z)),rL,'*','Color',ColorSeq(z,:),'LineWidth',1.5);
end
for z=1:length(LSet)
    PredDiff3=Rlower3(z);
    plot([1 N],PredDiff3*[1 1],'--','Color',.5*ColorSeq(z,:),'LineWidth',1.5);
end
xlim([1 N]); ylim([0 .6]);
set(gca,'XScale','log');
hl=legend(sprintf('$\\ell=$%.f',LSet(1)),sprintf('$\\ell=$%.f',LSet(2)),...
    sprintf('$\\ell=$%.f',LSet(3)), sprintf('$\\ell=$%.f',LSet(4)),...
    sprintf('$\\ell=$%.f',LSet(5)),'limits','theory prediction','inflection','location','Southeast');
set(hl,'interpreter','latex');
xlabel('Time','interpreter','latex');
ylabel('Mean growth rate','interpreter','latex');
title('Growth: Constant-to-fluctuating environments, $\beta=\beta_{Crit}/2$','interpreter','latex');


%% II. Adaptive l constant to fluctuating environment

adaptation='Variance'; adaptation='Proximity' %Toggle
rL=.05;                                 %rL=.5 %Toggle
rH=1;
p=.6;                                   p=.2 %TOGGLE
betaCrit=1+rH/rL
beta=betaCrit/2
pI=(beta-1)*rL/((beta-1)*rL+rH)
N=10^3;

method=NaN;
Hx=0;                                   Hx=1 %TOGGLE
NTraj=1;

lMin=3; %(lmin=2 to prevent trapping).
lMax=20; %TOGGLE lMax=20 and lMax=Inf

%LINEAR
lInf=lMax-(lMax-lMin)*abs(p-pI) %Longterm adaptive l
lconst=lMax-(lMax-lMin)*abs(Hx-pI)
LSet=[lMin round(lconst) min([lMax 20])];
%LSet=[5 10 15 25]

Hits=nan(1,length(LSet));

%Fixed memory:
StoreRAvgConst=cell(1,length(LSet));
%Variable memory:
StoreRAvgVariable=cell(1,length(LSet));
StoreLSeqAvgVariable=cell(1,length(LSet));
LSeqStorageVariable=nan(NTraj,N);

NIterate=10^3;
t=1;
for l=LSet
       %Constant
    RTotalsConst   =nan(NIterate,N+l);
    HitTotalsConst =nan(NIterate,N);
    LSeqTotalsConst=nan(NIterate,N);
        %Variable
    RTotalsVariable   =nan(NIterate,N+l);
    HitTotalsVariable =nan(NIterate,N);
    LSeqTotalsVariable=nan(NIterate,N);
    
    parfor z=1:NIterate
        [RConst HitConst RVariable HitVariable Env LSeqVariable] =...
            StochasticSensing1ConstantAndAdaptiveLPIProxv2(...
                rL, rH, beta, p, l, lMin, lMax, N, method,Hx,adaptation);    
        RTotalsConst(z,:)=RConst;
        HitTotalsConst(z,:)=HitConst;
        RTotalsVariable(z,:)=RVariable;
        HitTotalsVariable(z,:)=HitVariable;
        LSeqTotalsVariable(z,:)=LSeqVariable;
    end
    %Constant
    RAvgConst=mean(RTotalsConst,1);
    StoreRAvgConst{t}=RAvgConst;
    %Variable
    RAvgVariable=mean(RTotalsVariable,1);
    LSeqAvgVariable=mean(LSeqTotalsVariable,1);
    StoreRAvgVariable{t}=RAvgVariable;
    StoreLSeqAvgVariable{t}=LSeqAvgVariable;

    LSeqStorage{t}=LSeqTotalsVariable;
    t=t+1;
end

figure; hold on; box on;
for z=1:length(LSet)
    LSeqAvgVariable=StoreLSeqAvgVariable{z};
    plot(LSeqAvgVariable,'LineWidth',2);
end
plot([1:1:N], round(lInf*ones(1,N)),'k--','LineWidth',1.5);
LSeqTotalsVariable=LSeqStorage{2};
for z=1:NTraj
    plot(LSeqTotalsVariable(z,:),'LineWidth',2,'color',[ColorSeq(2,:) .25]);
end
set(gca,'XScale','log'); xlim([1 N]); ylim([1 lMax]);
hl=legend(sprintf('$N_0=N_\\textup{min}$'),sprintf('$N_0=N_c$'),...
    sprintf('$N_0=N_\\textup{Max}$'), sprintf('$N_\\infty$'),...
    sprintf('stochastic trajectories $N_0$=$N_c$'),'location','Southeast');
set(hl,'interpreter','latex');
xlabel('Time','interpreter','latex');
ylabel('Memory size $N_n$','interpreter','latex');
title(sprintf('%s-selected memory size $N$: constant-to-fluctuating environment ($p_I=$%.2f, $p_0=$%.2f, $\\beta=\\beta_{Crit}/2$, $\\ell_{Max}=$%.f)',adaptation,pI,p,lMax),'interpreter','latex');

%Comparative advantage of adaptive l versus fixed.
figure; hold on; box on;
ColorSeq=[0 .45 .74 ; .85 .33 .1; .93 .69 .13; .49 .18 .56; .47 .67 .19];
plot(0,0,'k-','linewidth',2);
plot(0,0,'k--','linewidth',2);
plot(0,0,'-','Color',ColorSeq(1,:),'linewidth',2);
plot(0,0,'-','Color',ColorSeq(2,:),'linewidth',2);
plot(0,0,'-','Color',ColorSeq(3,:),'linewidth',2);
plot(0,0,'k-.','Linewidth',1.5);
for z=1:3
    RAvgConst=StoreRAvgConst{z};
    plot([-LSet(z)+1:1:N],cumsum(RAvgConst)./[1:1:N+LSet(z)],'-','Color',ColorSeq(z,:),'LineWidth',2);
%     plot([-LSet(z)+1:1:N],RAvgConst,'-','Color',ColorSeq(z,:),'LineWidth',2);
end
for z=1:3
    RAvgVariable=StoreRAvgVariable{z};
    plot([-LSet(z)+1:1:N],cumsum(RAvgVariable)./[1:1:N+LSet(z)],'--','Color',.6*ColorSeq(z,:),'LineWidth',2);
%     plot([-LSet(z)+1:1:N],RAvgVariable,'--','Color',ColorSeq(z,:),'LineWidth',2);
end
set(gca,'XScale','log');
xlim([0 N]); %ylim([.25 .75]);
set(gca,'XScale','log');
hl=legend(sprintf('fixed $\\ell$'),sprintf('adaptive $\\ell$'),...
    sprintf('$\\ell_0=\\ell_\\textup{min}$'),sprintf('$\\ell_0=\\ell_c$'),...
    sprintf('$\\ell_0=\\ell_\\textup{Max}$'),'location','Southeast');
set(hl,'interpreter','latex');
xlabel('Time','interpreter','latex');
ylabel('Cumulative Averaged Growth','interpreter','latex');
title(sprintf('%s-selected adaptive vs. fixed memory growth: constant-to-fluctuating environment ($p_I=$%.2f, $p_0=$%.2f, $\\beta=\\beta_{Crit}/2$, $\\ell_{Max}=$%.f)',adaptation,pI,p,lMax),'interpreter','latex');


%% III. Adaptive l fluctuating to constant environment

adaptation='Variance'; adaptation='Proximity'
rL=.05;
rH=1;                                                       
p=.20; p=.5; p=.6                                   
betaCrit=1+rH/rL
beta=betaCrit/2 %Needs to be less than betaCrit
pI=(beta-1)*rL/((beta-1)*rL+rH);
N=1000;
method=ones(1,N);                               method=zeros(1,N) 
%Hx=1; %Defined later for this case
NTraj=10;

lMin=3;
lMax=20;

lInf=lMax-(lMax-lMin)*abs(method(1)-pI) %Longterm adaptive l
lconst=lMax-(lMax-lMin)*abs(p-pI) %adaptive l for previous state.

LSet=[lMin round(lconst) lMax]

Hits=nan(1,length(LSet));

%Fixed memory:
StoreRAvgConst=cell(1,length(LSet));
%Variable memory:
StoreRAvgVariable=cell(1,length(LSet));
StoreLSeqAvgVariable=cell(1,length(LSet));
LSeqStorageVariable=nan(NTraj,N);

NIterate=10^3;
t=1;
for l=LSet
        %Constant
    RTotalsConst   =nan(NIterate,N+l);
    HitTotalsConst =nan(NIterate,N);
    LSeqTotalsConst=nan(NIterate,N);
        %Variable
    RTotalsVariable   =nan(NIterate,N+l);
    HitTotalsVariable =nan(NIterate,N);
    LSeqTotalsVariable=nan(NIterate,N);
    
    for z=1:NIterate
        %Random
        Hx=rand(l,1)'<p;
        [RConst HitConst RVariable HitVariable Env LSeqVariable] =...
            StochasticSensing1ConstantAndAdaptiveLPIProxv2(...
                rL, rH, beta, p, l, lMin, lMax, N, method,Hx,adaptation);    
        RTotalsConst(z,:)=RConst;
        HitTotalsConst(z,:)=HitConst;
        RTotalsVariable(z,:)=RVariable;
        HitTotalsVariable(z,:)=HitVariable;
        LSeqTotalsVariable(z,:)=LSeqVariable;
    end
    %Constant
    RAvgConst=mean(RTotalsConst,1);
    StoreRAvgConst{t}=RAvgConst;
    %Variable
    RAvgVariable=mean(RTotalsVariable,1);
    LSeqAvgVariable=mean(LSeqTotalsVariable,1);
    StoreRAvgVariable{t}=RAvgVariable;
    StoreLSeqAvgVariable{t}=LSeqAvgVariable;

    LSeqStorage{t}=LSeqTotalsVariable;
    t=t+1;
end

figure; hold on; box on;
for z=1:length(LSet)
    LSeqAvgVariable=StoreLSeqAvgVariable{z};
    plot(LSeqAvgVariable,'LineWidth',2);
end
plot([1:1:N], round(lInf*ones(1,N)),'k--','LineWidth',1.5);
LSeqTotalsVariable=LSeqStorage{2};
for z=1:NTraj
    plot(LSeqTotalsVariable(z,:),'LineWidth',2,'color',[ColorSeq(2,:) .25]);
end
plot([1:1:N], round(lInf*ones(1,N)),'k--','LineWidth',1.5);

set(gca,'XScale','log'); xlim([1 N]); ylim([1 lMax]);
hl=legend(sprintf('$\\ell_0=\\ell_\\textup{min}$'),sprintf('$\\ell_0=\\ell_c$'),...
    sprintf('$\\ell_0=\\ell_\\textup{Max}$'), sprintf('$\\ell_\\infty$'),...
    sprintf('stochastic trajectories $\\ell_0$=$\\ell_c$'),'location','Southeast');
set(hl,'interpreter','latex');
xlabel('Time','interpreter','latex');
ylabel('Memory size $\ell$','interpreter','latex');
title(sprintf('%s-selected memory size $\\ell$: fluctuating-to-constant environment ($p_I=$%.2f, $p_0=$%.2f, $\\beta=\\beta_{Crit}/2$, $\\ell_{Max}=$%.f)',adaptation,pI,p,lMax),'interpreter','latex');

%Comparative advantage of adaptive l versus fixed.
figure; hold on; box on;
ColorSeq=[0 .45 .74 ; .85 .33 .1; .93 .69 .13; .49 .18 .56; .47 .67 .19];
plot(0,0,'k-','linewidth',2);
plot(0,0,'k--','linewidth',2);
plot(0,0,'-','Color',ColorSeq(1,:),'linewidth',2);
plot(0,0,'-','Color',ColorSeq(2,:),'linewidth',2);
plot(0,0,'-','Color',ColorSeq(3,:),'linewidth',2);
plot(0,0,'k-.','Linewidth',1.5);
for z=1:3
    RAvgConst=StoreRAvgConst{z};
   plot([-LSet(z)+1:1:N],cumsum(RAvgConst)./[1:1:N+LSet(z)],'-','Color',ColorSeq(z,:),'LineWidth',2);
end
for z=1:3
    RAvgVariable=StoreRAvgVariable{z};
    plot([-LSet(z)+1:1:N],cumsum(RAvgVariable)./[1:1:N+LSet(z)],'--','Color',.6*ColorSeq(z,:),'LineWidth',2);
end
set(gca,'XScale','log');
xlim([0 N]);
set(gca,'XScale','log');
hl=legend(sprintf('fixed $\\ell$'),sprintf('adaptive $\\ell$'),...
    sprintf('$\\ell_0=\\ell_\\textup{min}$'),sprintf('$\\ell_0=\\ell_c$'),...
    sprintf('$\\ell_0=\\ell_\\textup{Max}$'),'location','Southeast');
set(hl,'interpreter','latex');
xlabel('Time','interpreter','latex');
ylabel('Mean growth rate','interpreter','latex');
title(sprintf('%s-selected adaptive vs. fixed memory growth: fluctuating-to-constant environment ($p_I=$%.2f, $p_0=$%.2f, $\\beta=\\beta_{Crit}/2$, $\\ell_{Max}=$%.f)',adaptation,pI,p,lMax),'interpreter','latex');

%% Jensen's Inequality:

adaptation='Proximity'
rL=.1;  rL=0.01
rH=1;

betaCrit=1+rH/rL
beta=betaCrit/4
pI=(beta-1)*rL/((beta-1)*rL+rH)
N=10^3;
method=NaN;
NTraj=1;
lMin=3;
lMax=20; 

%I. Low constant:
p_CLow=0; p=p_CLow;
Hx_CLow=0; Hx=Hx_CLow;

lInf=lMax-(lMax-lMin)*abs(p-pI) %Longterm adaptive l
lconst=lMax-(lMax-lMin)*abs(Hx-pI)
LSet=[lMin round(lconst) min([lMax 20])]

Hits=nan(1,length(LSet));

%Fixed memory:
StoreRAvgConst=cell(1,length(LSet));
%Variable memory:
StoreRAvgVariable=cell(1,length(LSet));
StoreLSeqAvgVariable=cell(1,length(LSet));
LSeqStorageVariable=nan(NTraj,N);

NIterate=10^3;
t=1;
for l=LSet
        %Constant
    RTotalsConst   =nan(NIterate,N+l);
    HitTotalsConst =nan(NIterate,N);
    LSeqTotalsConst=nan(NIterate,N);
        %Variable
    RTotalsVariable   =nan(NIterate,N+l);
    HitTotalsVariable =nan(NIterate,N);
    LSeqTotalsVariable=nan(NIterate,N);
    
    parfor z=1:NIterate
        [RConst HitConst RVariable HitVariable Env LSeqVariable] =...
            StochasticSensing1ConstantAndAdaptiveLPIProxv2(...
                rL, rH, beta, p, l, lMin, lMax, N, method,Hx,adaptation);    
        RTotalsConst(z,:)=RConst;
        HitTotalsConst(z,:)=HitConst;
        RTotalsVariable(z,:)=RVariable;
        HitTotalsVariable(z,:)=HitVariable;
        LSeqTotalsVariable(z,:)=LSeqVariable;
    end
    %Constant
    RAvgConst=mean(RTotalsConst,1);
    StoreRAvgConst{t}=RAvgConst;
    %Variable
    RAvgVariable=mean(RTotalsVariable,1);
    LSeqAvgVariable=mean(LSeqTotalsVariable,1);
    StoreRAvgVariable{t}=RAvgVariable;
    StoreLSeqAvgVariable{t}=LSeqAvgVariable;

    LSeqStorage{t}=LSeqTotalsVariable;
    t=t+1;
end
RAvgConst_CLow=StoreRAvgConst{2}; RAvgConst_CLow=RAvgConst_CLow(LSet(2)+1:end);
RAvgVariable_CLow=StoreRAvgVariable{2}; RAvgVariable_CLow=RAvgVariable_CLow(LSet(2)+1:end);

%II. High constant:
p_CHigh=1; p=p_CHigh;
Hx_CHigh=1;Hx=Hx_CHigh;

lInf=lMax-(lMax-lMin)*abs(p-pI) %Longterm adaptive l
lconst=lMax-(lMax-lMin)*abs(Hx-pI)
LSet=[lMin round(lconst) min([lMax 20])]

Hits=nan(1,length(LSet));

%Fixed memory:
StoreRAvgConst=cell(1,length(LSet));
%Variable memory:
StoreRAvgVariable=cell(1,length(LSet));
StoreLSeqAvgVariable=cell(1,length(LSet));
LSeqStorageVariable=nan(NTraj,N);

NIterate=10^3;
t=1;
for l=LSet
    %Average Response
        %Constant
    RTotalsConst   =nan(NIterate,N+l);
    HitTotalsConst =nan(NIterate,N);
    LSeqTotalsConst=nan(NIterate,N);
        %Variable
    RTotalsVariable   =nan(NIterate,N+l);
    HitTotalsVariable =nan(NIterate,N);
    LSeqTotalsVariable=nan(NIterate,N);
    
    parfor z=1:NIterate
        [RConst HitConst RVariable HitVariable Env LSeqVariable] =...
            StochasticSensing1ConstantAndAdaptiveLPIProxv2(...
                rL, rH, beta, p, l, lMin, lMax, N, method,Hx,adaptation);    
        RTotalsConst(z,:)=RConst;
        HitTotalsConst(z,:)=HitConst;
        RTotalsVariable(z,:)=RVariable;
        HitTotalsVariable(z,:)=HitVariable;
        LSeqTotalsVariable(z,:)=LSeqVariable;
    end
    %Constant
    RAvgConst=mean(RTotalsConst,1);
    StoreRAvgConst{t}=RAvgConst;
    %Variable
    RAvgVariable=mean(RTotalsVariable,1);
    LSeqAvgVariable=mean(LSeqTotalsVariable,1);
    StoreRAvgVariable{t}=RAvgVariable;
    StoreLSeqAvgVariable{t}=LSeqAvgVariable;

    LSeqStorage{t}=LSeqTotalsVariable;
    t=t+1;
end
RAvgConst_CHigh=StoreRAvgConst{2}; RAvgConst_CHigh=RAvgConst_CHigh(LSet(2)+1:end);
RAvgVariable_CHigh=StoreRAvgVariable{2}; RAvgVariable_CHigh=RAvgVariable_CHigh(LSet(2)+1:end);

%IV. Fluctuating Avg
Delta=.1;
pLow=pI-Delta; pLow=pI;
pMed=pI; pMed=pI+Delta;
pHigh=pI+Delta; pHigh=pI+2*Delta;


P=[pI pI+Delta pI+2*Delta 1-Delta]; %P=[Delta pI pI+2*Delta 1-Delta];

for z=1:length(P)
    p=P(z)
    RAvgConst_CAverage(z,:) = RAvgConst_CLow + P(z)*(RAvgConst_CHigh-RAvgConst_CLow);
    RAvgVariable_CAverage(z,:) = RAvgVariable_CLow + P(z)*(RAvgVariable_CHigh-RAvgVariable_CLow);
    
    Hx=1; Hx=0; %TOGGLE
    
    t=1;
    for l=LSet
        %Average Response
        %Constant
        RTotalsConst   =nan(NIterate,N+l);
        HitTotalsConst =nan(NIterate,N);
        LSeqTotalsConst=nan(NIterate,N);
        %Variable
        RTotalsVariable   =nan(NIterate,N+l);
        HitTotalsVariable =nan(NIterate,N);
        LSeqTotalsVariable=nan(NIterate,N);
        
        parfor z=1:NIterate
            [RConst HitConst RVariable HitVariable Env LSeqVariable] =...
                StochasticSensing1ConstantAndAdaptiveLPIProxv2(...
                rL, rH, beta, p, l, lMin, lMax, N, method,Hx,adaptation);
            RTotalsConst(z,:)=RConst;
            HitTotalsConst(z,:)=HitConst;
            RTotalsVariable(z,:)=RVariable;
            HitTotalsVariable(z,:)=HitVariable;
            LSeqTotalsVariable(z,:)=LSeqVariable;
        end
        %Constant
        RAvgConst=mean(RTotalsConst,1);
        StoreRAvgConst{t}=RAvgConst;
        %Variable
        RAvgVariable=mean(RTotalsVariable,1);
        LSeqAvgVariable=mean(LSeqTotalsVariable,1);
        StoreRAvgVariable{t}=RAvgVariable;
        StoreLSeqAvgVariable{t}=LSeqAvgVariable;
        
        LSeqStorage{t}=LSeqTotalsVariable;
        t=t+1;
    end
    RAvgConst_Random=StoreRAvgConst{2}; RAvgConst_Random=RAvgConst_Random(LSet(2)+1:end);
    RAvgVariable_Random=StoreRAvgVariable{2}; RAvgVariable_Random=RAvgVariable_Random(LSet(2)+1:end);
    
    RAvgConst_RandomStore(z,:)=RAvgConst_Random;
    RAvgVariable_RandomStore(z,:)=RAvgVariable_Random;
end


%Comparative advantage of adaptive l versus fixed.
figure; hold on; box on;
ColorSeq=[0 .45 .74 ; .85 .33 .1; .93 .69 .13; .49 .18 .56; .47 .67 .19];
for z=1:length(P)
subplot(1,4,z);
if z==1 || z==3
plot(0,0,'k-','linewidth',2); hold on;
plot(0,0,'k--','linewidth',2);
plot(0,0,'-','Color',ColorSeq(1,:),'linewidth',2);
plot(0,0,'-','Color',ColorSeq(2,:),'linewidth',2);
plot(0,0,'-','Color',ColorSeq(3,:),'linewidth',2);
plot(0,0,'-','Color',ColorSeq(4,:),'linewidth',2);
end
plot([1:1:N],cumsum(RAvgConst_CLow)./[1:1:N],'-','Color',ColorSeq(1,:),'LineWidth',2); hold on;
plot([1:1:N],cumsum(RAvgConst_CHigh)./[1:1:N],'-','Color',ColorSeq(2,:),'LineWidth',2);
plot([1:1:N],cumsum(RAvgConst_CAverage(z,:))./[1:1:N],'-','Color',ColorSeq(3,:),'LineWidth',2);
plot([1:1:N],cumsum(RAvgConst_RandomStore(z,:))./[1:1:N],'-','Color',ColorSeq(4,:),'LineWidth',2);

plot([1:1:N],cumsum(RAvgVariable_CLow)./[1:1:N],'--','Color',.6*ColorSeq(1,:),'LineWidth',2);
plot([1:1:N],cumsum(RAvgVariable_CHigh)./[1:1:N],'--','Color',.6*ColorSeq(2,:),'LineWidth',2);
plot([1:1:N],cumsum(RAvgVariable_CAverage(z,:))./[1:1:N],'--','Color',.6*ColorSeq(3,:),'LineWidth',2);
plot([1:1:N],cumsum(RAvgVariable_RandomStore(z,:))./[1:1:N],'--','Color',.6*ColorSeq(4,:),'LineWidth',2);
set(gca,'XScale','log');
ylim([0 1]);
if z==1 || z==3
    hl=legend(sprintf('fixed $N$'),sprintf('adaptive $N$'),...
    sprintf('$C_\\textup{low}$'),sprintf('$C_\\textup{high}$'),...
    sprintf('$C_\\textup{avg}$'),sprintf('$p=C_\\textup{avg}$'),'location','East');
    set(hl,'interpreter','latex');
    ylabel('Cumulative Averaged Growth','interpreter','latex');
end
if z==1
    title(sprintf('$p=p_I=%.2f$',pI),'interpreter','latex');
elseif z==2
    title(sprintf('$p=%.2f$',pI+Delta),'interpreter','latex');
elseif z==3
    title(sprintf('$p=%.2f$',pI+2*Delta),'interpreter','latex');
    ylabel('Cumulative Averaged Growth','interpreter','latex');
elseif z==4
    title(sprintf('$p=%.2f$',1-Delta),'interpreter','latex');
end
xlabel('Time','interpreter','latex');
set(gca,'FontSize',16);
end
h=sgtitle(sprintf('Cumulative averaged growth: \n Constant and rapidly fluctuating nutrient environments'))
h.Interpreter='latex';
h.FontSize=20

%% Jensen's across 2-parameter perturbation:
adaptation='Proximity';
rH=1; rL=.01;
betaCrit=1+rH/rL
method=NaN;
%Range beta: 1 to betaCrit
%Range p=0.01 to 0.99:
%plot Jensen Difference:
BETA=1:floor(betaCrit);
P=0.01:0.01:0.99;
NIterate=10^2;
NTraj=1;
N=10^3;
tic;
StorageRAvgConstLow=nan(length(P),length(BETA));
StorageRAvgVariableLow=nan(length(P),length(BETA));

StorageRAvgConstHigh=nan(length(P),length(BETA));
StorageRAvgVariableHigh=nan(length(P),length(BETA));

StorageRAvgConstAverage=nan(length(P),length(BETA));
StorageRAvgVariableAverage=nan(length(P),length(BETA));

StorageRAvgConstRandom=nan(length(P),length(BETA));
StorageRAvgVariableRandom=nan(length(P),length(BETA));
lMin=2; lMax=20;
for z1=1:length(P)
    toc; z1
    for z2=1:length(BETA)
        p=P(z1);
        beta=BETA(z2);
        pI=(beta-1)*rL/((beta-1)*rL+rH);
        %I. Low
        p_CLow=0; p=p_CLow;
        Hx_CLow=0; Hx=Hx_CLow;
        %lInf=lMax-(lMax-lMin)*abs(p-pI) %Longterm adaptive l
        lconst=lMax-(lMax-lMin)*abs(Hx-pI);
        LSet=round(lconst); %LSet=[lMin round(lconst) min([lMax 20])];
        Hits=nan(1,length(LSet));
        %Fixed memory:
        StoreRAvgConst=cell(1,length(LSet));
        %Variable memory:
        StoreRAvgVariable=cell(1,length(LSet));
        StoreLSeqAvgVariable=cell(1,length(LSet));
        LSeqStorageVariable=nan(NTraj,N);
        t=1;
        for l=LSet
            %Average Response
            %Constant
            RTotalsConst   =nan(NIterate,N+l);
            HitTotalsConst =nan(NIterate,N);
            LSeqTotalsConst=nan(NIterate,N);
            %Variable
            RTotalsVariable   =nan(NIterate,N+l);
            HitTotalsVariable =nan(NIterate,N);
            LSeqTotalsVariable=nan(NIterate,N);
            parfor z=1:NIterate
                [RConst HitConst RVariable HitVariable Env LSeqVariable] =...
                    StochasticSensing1ConstantAndAdaptiveLPIProxv2(...
                    rL, rH, beta, p, l, lMin, lMax, N, method,Hx,adaptation);
                RTotalsConst(z,:)=RConst;
                HitTotalsConst(z,:)=HitConst;
                RTotalsVariable(z,:)=RVariable;
                HitTotalsVariable(z,:)=HitVariable;
                LSeqTotalsVariable(z,:)=LSeqVariable;
            end
            %Constant
            RAvgConst=mean(RTotalsConst,1);
            StoreRAvgConst{t}=RAvgConst;
            %Variable
            RAvgVariable=mean(RTotalsVariable,1);
            LSeqAvgVariable=mean(LSeqTotalsVariable,1);
            StoreRAvgVariable{t}=RAvgVariable;
            StoreLSeqAvgVariable{t}=LSeqAvgVariable;
            LSeqStorage{t}=LSeqTotalsVariable;
            t=t+1;
        end
        RAvgConst_CLow=StoreRAvgConst{1}; RAvgConst_CLow=RAvgConst_CLow(LSet(1)+1:end);
        RAvgVariable_CLow=StoreRAvgVariable{1}; RAvgVariable_CLow=RAvgVariable_CLow(LSet(1)+1:end);
        StorageRAvgConstLow(z1,z2)=RAvgConst_CLow(end);
        StorageRAvgVariableLow(z1,z2)=RAvgVariable_CLow(end);

        %II. High
        p_CHigh=1; p=p_CHigh;
        Hx_CHigh=1;Hx=Hx_CHigh;
        %lInf=lMax-(lMax-lMin)*abs(p-pI) %Longterm adaptive l
        lconst=lMax-(lMax-lMin)*abs(Hx-pI);
        LSet=round(lconst); %LSet=[lMin round(lconst) min([lMax 20])]
        Hits=nan(1,length(LSet));
        %Fixed memory:
        StoreRAvgConst=cell(1,length(LSet));
        %Variable memory:
        StoreRAvgVariable=cell(1,length(LSet));
        StoreLSeqAvgVariable=cell(1,length(LSet));
        LSeqStorageVariable=nan(NTraj,N);
        t=1;
        for l=LSet
            %Average Response
            %Constant
            RTotalsConst   =nan(NIterate,N+l);
            HitTotalsConst =nan(NIterate,N);
            LSeqTotalsConst=nan(NIterate,N);
            %Variable
            RTotalsVariable   =nan(NIterate,N+l);
            HitTotalsVariable =nan(NIterate,N);
            LSeqTotalsVariable=nan(NIterate,N);
            parfor z=1:NIterate
                [RConst HitConst RVariable HitVariable Env LSeqVariable] =...
                    StochasticSensing1ConstantAndAdaptiveLPIProxv2(...
                    rL, rH, beta, p, l, lMin, lMax, N, method,Hx,adaptation);
                RTotalsConst(z,:)=RConst;
                HitTotalsConst(z,:)=HitConst;
                RTotalsVariable(z,:)=RVariable;
                HitTotalsVariable(z,:)=HitVariable;
                LSeqTotalsVariable(z,:)=LSeqVariable;
            end
            %Constant
            RAvgConst=mean(RTotalsConst,1);
            StoreRAvgConst{t}=RAvgConst;
            %Variable
            RAvgVariable=mean(RTotalsVariable,1);
            LSeqAvgVariable=mean(LSeqTotalsVariable,1);
            StoreRAvgVariable{t}=RAvgVariable;
            StoreLSeqAvgVariable{t}=LSeqAvgVariable;
            LSeqStorage{t}=LSeqTotalsVariable;
            t=t+1;
        end
        RAvgConst_CHigh=StoreRAvgConst{1}; RAvgConst_CHigh=RAvgConst_CHigh(LSet(1)+1:end);
        RAvgVariable_CHigh=StoreRAvgVariable{1}; RAvgVariable_CHigh=RAvgVariable_CHigh(LSet(1)+1:end);
        StorageRAvgConstHigh(z1,z2)=RAvgConst_CHigh(end);
        StorageRAvgVariableHigh(z1,z2)=RAvgVariable_CHigh(end);


        StorageRAvgConstAverage(z1,z2)=StorageRAvgConstLow(z1,z2)+P(z1)*(StorageRAvgConstHigh(z1,z2)-StorageRAvgConstLow(z1,z2));
        StorageRAvgVariableAverage(z1,z2)=StorageRAvgVariableLow(z1,z2)+P(z1)*(StorageRAvgVariableHigh(z1,z2)-StorageRAvgVariableLow(z1,z2));
        %III. Fluctuations
        p=P(z1);
        beta=BETA(z2);
        
        Hx=1; Hx=0; %TOGGLE
        lconst=lMax-(lMax-lMin)*abs(Hx-pI);
        LSet=round(lconst); %LSet=[lMin round(lconst) min([lMax 20])]
        t=1;
        for l=LSet
            %Average Response
            %Constant
            RTotalsConst   =nan(NIterate,N+l);
            HitTotalsConst =nan(NIterate,N);
            LSeqTotalsConst=nan(NIterate,N);
            %Variable
            RTotalsVariable   =nan(NIterate,N+l);
            HitTotalsVariable =nan(NIterate,N);
            LSeqTotalsVariable=nan(NIterate,N);
            parfor z=1:NIterate
                [RConst HitConst RVariable HitVariable Env LSeqVariable] =...
                    StochasticSensing1ConstantAndAdaptiveLPIProxv2(...
                    rL, rH, beta, p, l, lMin, lMax, N, method,Hx,adaptation);
                RTotalsConst(z,:)=RConst;
                HitTotalsConst(z,:)=HitConst;
                RTotalsVariable(z,:)=RVariable;
                HitTotalsVariable(z,:)=HitVariable;
                LSeqTotalsVariable(z,:)=LSeqVariable;
            end
            %Constant
            RAvgConst=mean(RTotalsConst,1);
            StoreRAvgConst{t}=RAvgConst;
            %Variable
            RAvgVariable=mean(RTotalsVariable,1);
            LSeqAvgVariable=mean(LSeqTotalsVariable,1);
            StoreRAvgVariable{t}=RAvgVariable;
            StoreLSeqAvgVariable{t}=LSeqAvgVariable;
            LSeqStorage{t}=LSeqTotalsVariable;
            t=t+1;
        end
        RAvgConst_Random=StoreRAvgConst{1}; RAvgConst_Random=RAvgConst_Random(LSet(1)+1:end);
        RAvgVariable_Random=StoreRAvgVariable{1}; RAvgVariable_Random=RAvgVariable_Random(LSet(1)+1:end);
        StorageRAvgConstRandom(z1,z2)=RAvgConst_Random(end);
        StorageRAvgVariableRandom(z1,z2)=RAvgVariable_Random(end);
    end
end

addpath /Users/jasontgeorge/Documents/MATLAB/Functions;

[X Y]=meshgrid(BETA/betaCrit,P);

%Simulated stochastic
figure; hold on; box on; grid on;
surfc(X,Y,StorageRAvgVariableRandom,'EdgeAlpha',0);
xlabel(sprintf('Relative growth coefficient $\\beta/\\beta_\\textup{c}$'),'interpreter','latex');
ylabel(sprintf('Convex combination parameter $\\varphi$'),'interpreter','latex');
zlabel(sprintf('$R_\\textup{rand}$'),'interpreter','latex');
xlim([0 1]); ylim([0 1]); zlim([0 1]);
title(sprintf('Equilibrium growth rate in a random environment $R_\\textup{rand}$'),'interpreter','latex');
set(gca,'FontSize',14);

%Averaged
figure; hold on; box on; grid on;
surfc(X,Y,StorageRAvgVariableAverage,'EdgeAlpha',0);
colorbar;
xlabel(sprintf('Relative growth coefficient $\\beta/\\beta_\\textup{c}$'),'interpreter','latex');
ylabel(sprintf('Convex combination parameter $\\varphi$'),'interpreter','latex');
zlabel(sprintf('$R_\\textup{avg}$'),'interpreter','latex');
xlim([0 1]); ylim([0 1]); zlim([0 1]);
title(sprintf('Averaged equilibrium growth rate of constant environments $R_\\textup{avg}$'),'interpreter','latex');
set(gca,'FontSize',14);

BETA;
pI=(BETA-1)*rL./((BETA-1)*rL+rH);
MaxDiff=(BETA-1)*rL*rH./((BETA-1)*rL+rH);

figure; hold on; box on; grid on;
plot3(0,0,0,'r','LineWidth',8);
surfc(X,Y,StorageRAvgVariableAverage - StorageRAvgVariableRandom,'FaceAlpha',0.75,'EdgeAlpha',.5);
s.EdgeColor = 'none';
plot3(BETA/betaCrit,pI,MaxDiff,'r','LineWidth',8);
colorbar;
view([-1 -1 1]);
xlabel(sprintf('$\\beta/\\beta_\\textup{c}$'),'interpreter','latex');
ylabel(sprintf('$\\varphi$'),'interpreter','latex');
zlabel(sprintf('$D_\\textup{tot}=R_\\textup{avg} -  R_\\textup{rand}$'),'interpreter','latex');
xlim([0 1]); ylim([0 1]); %zlim([0 1]);
title(sprintf('$D_\\textup{tot}$ vs.  $\\varphi$ and $\\beta / \\beta_\\textup{c}$'),'interpreter','latex');
set(gca,'FontSize',20)
leg1=legend(sprintf('Predicted maximal \n intrinsic deficit \n $D_I(p_I)$'),'Location','West')
set(leg1,'Interpreter','latex');
axis equal

figure; hold on; box on; grid on;
surfc(X,Y,StorageRAvgVariableAverage - StorageRAvgVariableRandom);
colorbar;
xlabel(sprintf('Relative growth efficiency coefficient $\\beta/\\beta_\\textup{c}$'),'interpreter','latex');
ylabel(sprintf('Convex combination parameter $\\varphi$'),'interpreter','latex');
xlim([0 1]); ylim([0 1]); %zlim([0 .5]);
title(sprintf('$R_\\textup{avg} -  R_\\textup{rand}$ vs.  $\\varphi$ and $\\beta / \\beta_\\textup{c}$'),'interpreter','latex');
set(gca,'FontSize',16)


%% 10/30/21: Same as above, but limit memory M=3 to show non-linear exogenous deficiency
% Jensen's across 2-parameter perturbation:
adaptation='Proximity';
rH=1; rL=.01;
betaCrit=1+rH/rL
method=NaN;
BETA=1:floor(betaCrit);
P=0.01:0.01:0.99;
NIterate=10^3;
NTraj=1;
N=              10^3;
tic;
StorageRAvgConstLow=nan(length(P),length(BETA));
StorageRAvgVariableLow=nan(length(P),length(BETA));

StorageRAvgConstHigh=nan(length(P),length(BETA));
StorageRAvgVariableHigh=nan(length(P),length(BETA));

StorageRAvgConstAverage=nan(length(P),length(BETA));
StorageRAvgVariableAverage=nan(length(P),length(BETA));

StorageRAvgConstRandom=nan(length(P),length(BETA));
StorageRAvgVariableRandom=nan(length(P),length(BETA));
lMin=3; lMax=3;                                                         %lMin=11; lMax=11;
for z1=1:length(P)
    toc; z1
    for z2=1:length(BETA)
        p=P(z1);
        beta=BETA(z2);
        pI=(beta-1)*rL/((beta-1)*rL+rH);
        %I. Low
        p_CLow=0; p=p_CLow;
        Hx_CLow=0; Hx=Hx_CLow;
        %lInf=lMax-(lMax-lMin)*abs(p-pI) %Longterm adaptive l
        lconst=lMax-(lMax-lMin)*abs(Hx-pI);
        LSet=round(lconst); %LSet=[lMin round(lconst) min([lMax 20])];
        Hits=nan(1,length(LSet));
        %Fixed memory:
        StoreRAvgConst=cell(1,length(LSet));
        %Variable memory:
        StoreRAvgVariable=cell(1,length(LSet));
        StoreLSeqAvgVariable=cell(1,length(LSet));
        LSeqStorageVariable=nan(NTraj,N);
        t=1;
        for l=LSet
            %Average Response
            %Constant
            RTotalsConst   =nan(NIterate,N+l);
            HitTotalsConst =nan(NIterate,N);
            LSeqTotalsConst=nan(NIterate,N);
            %Variable
            RTotalsVariable   =nan(NIterate,N+l);
            HitTotalsVariable =nan(NIterate,N);
            LSeqTotalsVariable=nan(NIterate,N);
            parfor z=1:NIterate
                [RConst HitConst RVariable HitVariable Env LSeqVariable] =...
                    StochasticSensing1ConstantAndAdaptiveLPIProxv2(...
                    rL, rH, beta, p, l, lMin, lMax, N, method,Hx,adaptation);
                RTotalsConst(z,:)=RConst;
                HitTotalsConst(z,:)=HitConst;
                RTotalsVariable(z,:)=RVariable;
                HitTotalsVariable(z,:)=HitVariable;
                LSeqTotalsVariable(z,:)=LSeqVariable;
            end
            %Constant
            RAvgConst=mean(RTotalsConst,1);
            StoreRAvgConst{t}=RAvgConst;
            %Variable
            RAvgVariable=mean(RTotalsVariable,1);
            LSeqAvgVariable=mean(LSeqTotalsVariable,1);
            StoreRAvgVariable{t}=RAvgVariable;
            StoreLSeqAvgVariable{t}=LSeqAvgVariable;
            LSeqStorage{t}=LSeqTotalsVariable;
            t=t+1;
        end
        RAvgConst_CLow=StoreRAvgConst{1}; RAvgConst_CLow=RAvgConst_CLow(LSet(1)+1:end);
        RAvgVariable_CLow=StoreRAvgVariable{1}; RAvgVariable_CLow=RAvgVariable_CLow(LSet(1)+1:end);
        StorageRAvgConstLow(z1,z2)=RAvgConst_CLow(end);
        StorageRAvgVariableLow(z1,z2)=RAvgVariable_CLow(end);

        %II. High
        p_CHigh=1; p=p_CHigh;
        Hx_CHigh=1;Hx=Hx_CHigh;
        %lInf=lMax-(lMax-lMin)*abs(p-pI) %Longterm adaptive l
        lconst=lMax-(lMax-lMin)*abs(Hx-pI);
        LSet=round(lconst); %LSet=[lMin round(lconst) min([lMax 20])]
        Hits=nan(1,length(LSet));
        %Fixed memory:
        StoreRAvgConst=cell(1,length(LSet));
        %Variable memory:
        StoreRAvgVariable=cell(1,length(LSet));
        StoreLSeqAvgVariable=cell(1,length(LSet));
        LSeqStorageVariable=nan(NTraj,N);
        t=1;
        for l=LSet
            %Constant
            RTotalsConst   =nan(NIterate,N+l);
            HitTotalsConst =nan(NIterate,N);
            LSeqTotalsConst=nan(NIterate,N);
            %Variable
            RTotalsVariable   =nan(NIterate,N+l);
            HitTotalsVariable =nan(NIterate,N);
            LSeqTotalsVariable=nan(NIterate,N);
            parfor z=1:NIterate
                [RConst HitConst RVariable HitVariable Env LSeqVariable] =...
                    StochasticSensing1ConstantAndAdaptiveLPIProxv2(...
                    rL, rH, beta, p, l, lMin, lMax, N, method,Hx,adaptation);
                RTotalsConst(z,:)=RConst;
                HitTotalsConst(z,:)=HitConst;
                RTotalsVariable(z,:)=RVariable;
                HitTotalsVariable(z,:)=HitVariable;
                LSeqTotalsVariable(z,:)=LSeqVariable;
            end
            %Constant
            RAvgConst=mean(RTotalsConst,1);
            StoreRAvgConst{t}=RAvgConst;
            %Variable
            RAvgVariable=mean(RTotalsVariable,1);
            LSeqAvgVariable=mean(LSeqTotalsVariable,1);
            StoreRAvgVariable{t}=RAvgVariable;
            StoreLSeqAvgVariable{t}=LSeqAvgVariable;
            LSeqStorage{t}=LSeqTotalsVariable;
            t=t+1;
        end
        RAvgConst_CHigh=StoreRAvgConst{1}; RAvgConst_CHigh=RAvgConst_CHigh(LSet(1)+1:end);
        RAvgVariable_CHigh=StoreRAvgVariable{1}; RAvgVariable_CHigh=RAvgVariable_CHigh(LSet(1)+1:end);
        StorageRAvgConstHigh(z1,z2)=RAvgConst_CHigh(end);
        StorageRAvgVariableHigh(z1,z2)=RAvgVariable_CHigh(end);


        StorageRAvgConstAverage(z1,z2)=StorageRAvgConstLow(z1,z2)+P(z1)*(StorageRAvgConstHigh(z1,z2)-StorageRAvgConstLow(z1,z2));
        StorageRAvgVariableAverage(z1,z2)=StorageRAvgVariableLow(z1,z2)+P(z1)*(StorageRAvgVariableHigh(z1,z2)-StorageRAvgVariableLow(z1,z2));
        %III. Fluctuations
        p=P(z1);
        beta=BETA(z2);
        
        Hx=1; Hx=0; %TOGGLE
        lconst=lMax-(lMax-lMin)*abs(Hx-pI);
        LSet=round(lconst); %LSet=[lMin round(lconst) min([lMax 20])]
        t=1;
        for l=LSet
            %Average Response
            %Constant
            RTotalsConst   =nan(NIterate,N+l);
            HitTotalsConst =nan(NIterate,N);
            LSeqTotalsConst=nan(NIterate,N);
            %Variable
            RTotalsVariable   =nan(NIterate,N+l);
            HitTotalsVariable =nan(NIterate,N);
            LSeqTotalsVariable=nan(NIterate,N);
            parfor z=1:NIterate
                [RConst HitConst RVariable HitVariable Env LSeqVariable] =...
                    StochasticSensing1ConstantAndAdaptiveLPIProxv2(...
                    rL, rH, beta, p, l, lMin, lMax, N, method,Hx,adaptation);
                RTotalsConst(z,:)=RConst;
                HitTotalsConst(z,:)=HitConst;
                RTotalsVariable(z,:)=RVariable;
                HitTotalsVariable(z,:)=HitVariable;
                LSeqTotalsVariable(z,:)=LSeqVariable;
            end
            %Constant
            RAvgConst=mean(RTotalsConst,1);
            StoreRAvgConst{t}=RAvgConst;
            %Variable
            RAvgVariable=mean(RTotalsVariable,1);
            LSeqAvgVariable=mean(LSeqTotalsVariable,1);
            StoreRAvgVariable{t}=RAvgVariable;
            StoreLSeqAvgVariable{t}=LSeqAvgVariable;
            LSeqStorage{t}=LSeqTotalsVariable;
            t=t+1;
        end
        RAvgConst_Random=StoreRAvgConst{1}; RAvgConst_Random=RAvgConst_Random(LSet(1)+1:end);
        RAvgVariable_Random=StoreRAvgVariable{1}; RAvgVariable_Random=RAvgVariable_Random(LSet(1)+1:end);
        StorageRAvgConstRandom(z1,z2)=RAvgConst_Random(end);
        StorageRAvgVariableRandom(z1,z2)=RAvgVariable_Random(end);
    end
end

[X Y]=meshgrid(BETA/betaCrit,P);


BETA;
pI=(BETA-1)*rL./((BETA-1)*rL+rH);
MaxDiff=(BETA-1)*rL*rH./((BETA-1)*rL+rH);

%Numerically calculate max:
xMaxLow=nan(length(BETA),1);
xMaxHigh=nan(length(BETA),1);
zMaxLow=nan(length(BETA),1);
zMaxHigh=nan(length(BETA),1);
for j=1:length(BETA)%for each Beta
    %consider Diff(\phi), and maximize over \phi
    M=2; %M=11 %Try out larger M case!
    %Maximize: DE * normcdf((pI-x)*sqrt(M/(x*(1-x)),0,1) * (x*rH - (1-x)*(BETA(j)-1)*rL)
    %To be maximized on x\in[0,pI]
    DILow = @(x)  -(rH*x + ...
        abs((1-x).*(BETA(j)-1)*rL-x.*rH) .* (1 - normcdf((pI(j)-x).*sqrt(M./(x.*(1-x))),0,1) ) );
    %To be maximized on x\in[pI,1]
    DIHigh = @(x) -((BETA(j)-1)*rL.*(1-x) + ...
        abs((1-x).*(BETA(j)-1)*rL-x.*rH) .* (normcdf((pI(j)-x).*sqrt(M./(x.*(1-x))),0,1) ) );

    xMaxLow(j)  = fminbnd(DILow,0,pI(j));
    xMaxHigh(j) = fminbnd(DIHigh,pI(j),1);
    
    zMaxLow(j)  = -DILow(xMaxLow(j));
    zMaxHigh(j) = -DIHigh(xMaxHigh(j));
end
j=101;
DILow = @(x) -(rH*x + abs((1-x).*(BETA(j)-1)*rL-x.*rH) .* (1 - normcdf((pI(j)-x).*sqrt(M./(x.*(1-x))),0,1) ) );
DIHigh = @(x) -((BETA(j)-1)*rL.*(1-x) + abs((1-x).*(BETA(j)-1)*rL-x.*rH) .* (normcdf((pI(j)-x).*sqrt(M./(x.*(1-x))),0,1) ) ); 
DIEnd =-[DILow([0:.01:pI(j)]) DIHigh([pI(j):.01:1])]
j=50;
DILow = @(x) -(rH*x + abs((1-x).*(BETA(j)-1)*rL-x.*rH) .* (1 - normcdf((pI(j)-x).*sqrt(M./(x.*(1-x))),0,1) ) );
DIHigh = @(x) -((BETA(j)-1)*rL.*(1-x) + abs((1-x).*(BETA(j)-1)*rL-x.*rH) .* (normcdf((pI(j)-x).*sqrt(M./(x.*(1-x))),0,1) ) ); 
DIMid =-[DILow([0:.01:pI(j)]) DIHigh([pI(j):.01:1])] 
j=10;
DILow = @(x) -(rH*x + abs((1-x).*(BETA(j)-1)*rL-x.*rH) .* (1 - normcdf((pI(j)-x).*sqrt(M./(x.*(1-x))),0,1) ) );
DIHigh = @(x) -((BETA(j)-1)*rL.*(1-x) + abs((1-x).*(BETA(j)-1)*rL-x.*rH) .* (normcdf((pI(j)-x).*sqrt(M./(x.*(1-x))),0,1) ) ); 
DIEarly =-[DILow([0:.01:pI(j)]) DIHigh([pI(j):.01:1])] 


%TROUBLE-SHOOTING SLICES 

figure; hold on; box on; grid on;
plot3(-1,-1,-1,'LineWidth',2,'Color',[0.72 0.27 1.00]);
plot3(-1,-1,-1,'k-','LineWidth',4);
plot3(-1,-1,-1,'-','LineWidth',8,'Color',[1 0 0 1]);
plot3(-1,-1,-1,'-','LineWidth',8,'Color',[0.00,0.45,0.74 1]);
%Beta=BetaCrit so that pI=.5;
j=101; 
plot3(BETA(j)/betaCrit*ones(length(DIEnd),1), [0:.01:pI(j) pI(j):.01:1], DIEnd,'k-','LineWidth',4);
plot3(BETA(j)/betaCrit*ones(length(P),1),P,StorageRAvgVariableAverage(:,j) - StorageRAvgVariableRandom(:,j),'-','LineWidth',2,'Color',[0.72 0.27 1.00]);
%Beta=BetaCrit/2 so that pI=.5;
j=50;
plot3(BETA(j)/betaCrit*ones(length(DIMid),1), [0:.01:pI(j) pI(j):.01:1],DIMid,'k-','LineWidth',4);
plot3(BETA(j)/betaCrit*ones(length(P),1),P,StorageRAvgVariableAverage(:,j) - StorageRAvgVariableRandom(:,j),'-','LineWidth',2,'Color',[0.72 0.27 1.00]);
j=10;
plot3(BETA(j)/betaCrit*ones(length(DIEarly),1), [0:.01:pI(j) pI(j):.01:1],DIEarly,'k-','LineWidth',4);
plot3(BETA(j)/betaCrit*ones(length(P),1),P,StorageRAvgVariableAverage(:,j) - StorageRAvgVariableRandom(:,j),'-','LineWidth',2,'Color',[0.72 0.27 1.00]);
view([-45 35])
%Plot DI at PI
plot3(BETA/betaCrit,pI,MaxDiff,'-','LineWidth',8,'Color',[1 0 0 1]);
%Shadow
plot3(BETA/betaCrit,pI,zeros(length(MaxDiff),1),'-','LineWidth',8,'Color',[1 0 0 .2]);
plot3(BETA/betaCrit,ones(length(pI),1),MaxDiff,'-','LineWidth',8,'Color',[1 0 0 .2]);
%Plot Predicted Max
plot3(BETA/betaCrit,xMaxHigh,zMaxHigh,'LineWidth',8,'Color',[0.00,0.45,0.74 1]);
%Shadow
plot3(BETA/betaCrit,xMaxHigh,zeros(length(zMaxHigh),1),'r','LineWidth',8,'Color',[0.00,0.45,0.74 .2]);
plot3(BETA/betaCrit,ones(length(xMaxHigh),1),zMaxHigh,'r','LineWidth',8,'Color',[0.00,0.45,0.74 .2]);
xlim([0 1]); ylim([0 1]);
zlim([-.1 .6]);
xlabel(sprintf('$\\beta/\\beta_\\textup{c}$'),'interpreter','latex');
ylabel(sprintf('$p$'),'interpreter','latex');
zlabel(sprintf('Expected $D_\\textup{tot}$'),'interpreter','latex');
title(sprintf('Expected $D_\\textup{tot}$ vs.  $p$ and $\\beta / \\beta_\\textup{c}$'),'interpreter','latex');
set(gca,'FontSize',20)
leg=legend('Simulations','Theory',sprintf('$D_I(p_I)$'),sprintf('$\\max D_\\textup{tot}(p,\\beta)$'),'Location','West');
set(leg,'Interpreter','latex');

