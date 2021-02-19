%% Demonstration program for flexible project management
%Author: Zsolt T. Kosztyan Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Demonstration program for flexible project management
%Please cite as our works as:

%% Initialisation

close all
clear all

%% Genarate PDM

N=30      % Number of tasks
ff=0.3    % Flexiblity factor (between 0 to 1)
cf=3      % Connectivity factor (0,1,...) ~ number of parallel paths
mTD=20    % Max value of durations (must be positive)
mCD=30    % Max value of cost demands (must be positive)
mRD=20    % Max value of resource demands (must be positive)
w=2       % Number of completion mode
nR=2      % Number of resources
nW=0      % Number of unplanned tasks
scale=1.6 % Scale parameter

is_quality=0 %0 no quality domain, 1 use quality domain

if is_quality==0
    PDM=generatepdm_scale(N,ff,cf,mTD,mCD,mRD,w,nR,nW,scale);
else
    PDM=generatepdmq_scale(N,ff,cf,mTD,mCD,mRD,w,nR,nW,scale);
end

%% Draw topology

figure('Name','Topology of networks','units','normalized','outerposition',[0 0 1 1])
subplot(1,3,1)
LD=triu(PDM(:,1:N),1);
ld=LD;
ld(ld<1)=0;
G=digraph(LD);
g=digraph(ld);
h=plot(G);
Tasks=1:N;
highlight(h,g,'EdgeColor','r','LineWidth',1.5);
highlight(h,floor(Tasks(diag(PDM(:,1:N))==1)),'NodeColor','r')
title('Labelled graph');
subplot(1,3,2)
plot(G);
title('All tasks and dependencies');

subplot(1,3,3)
plot(digraph(ld(diag(PDM(:,1:N))==1,diag(PDM(:,1:N))==1)),'r');
title('Mandatory tasks and strict dependencies');


%% Desriptive statistics

figure('Name','Descriptive statistics of demands','units','normalized','outerposition',[0 0 1 1])

TPTmin=percentt(PDM,w,0);
TPTmax=percentt(PDM,w,1);
TPCmin=percentc(PDM,w,0);
TPCmax=percentc(PDM,w,1);
TPSmin=percents(PDM,0);
TPSmax=percents(PDM,1);

if is_quality==0
    subplot(1,4,1)
    bar(categorical({'TPT_{min}','TPT_{max}'}),[TPTmin,TPTmax])
    title('Duration')
    
    subplot(1,4,2)
    bar(categorical({'TPC_{min}','TPC_{max}'}),[TPCmin,TPCmax])    
    title('Cost demands')
    
    subplot(1,4,3)
    c=categorical({});
    TPRmin=percentr(PDM,w,0);
    TPRmax=percentr(PDM,w,1);
    TPR=[TPRmin,TPRmax];
    for j=1:size(TPR,2)
         c(j,1)=strcat('TPR_',num2str(j),'_{min}');
         c(j,2)=strcat('TPR_',num2str(j),'_{max}');
    end     
    bar(c,TPR)
    title('Resource demands')
    
    subplot(1,4,4)
    bar(categorical({'TPS_{min}','TPS_{max}'}),[TPSmin,TPSmax])    
    title('Score/scope demands')

else
    
    subplot(1,5,1)
    bar(categorical({'TPT_{min}','TPT_{max}'}),[TPTmin,TPTmax])
    title('Duration')
    
    subplot(1,5,2)
    bar(categorical({'TPC_{min}','TPC_{max}'}),[TPCmin,TPCmax])    
    title('Cost demands')
    
    subplot(1,5,3)
    TPQmin=percentq(PDM,w,0);
    TPQmax=percentq(PDM,w,1);

    bar(categorical({'TPQ_{min}','TPQ_{max}'}),[TPQmin,TPQmax])    
    title('Quality parameters')
    
    subplot(1,5,4)
    c=categorical({});
    TPRmin=percentqr(PDM,w,0);
    TPRmax=percentqr(PDM,w,1);
    TPR=[TPRmin,TPRmax];
    for j=1:size(TPR,2)
         c(j,1)=strcat('TPR_',num2str(j),'_{min}');
         c(j,2)=strcat('TPR_',num2str(j),'_{max}');
    end    
    bar(c,TPR)
    title('Resource demands')

    subplot(1,5,5)
    bar(categorical({'TPS_{min}','TPS_{max}'}),[TPSmin,TPSmax])    
    title('Score/scope demands')
end

%% Optimization

typefcn=1 % {0=maxTPQ,} 1=minTPT, 2=minTPC, 3=maxTPS,{4=minUF,} ~ composite

pT=0.85;
pC=0.85;
pQ=0.30;
pS=0.30;
pR=1.00;

Ct=percentt(PDM,w,pT);
Cc=percentc(PDM,w,pT);
Cs=percents(PDM,pS);

if is_quality==0
    Cr=percentr(PDM,w,pR)';
    CONS=[Ct,Cc,Cr,Cs]
    Select=10 %RC-(H)CTCTP
else
    Cq=percentq(PDM,w,pC);
    Cr=percentqr(PDM,w,pR)';
    CONS=[Ct,Cc,Cq,Cr,Cs]
    Select=17 %RC-(H)CTQCTP
end

PSM_tpma=mtpma(PDM,CONS,Select,typefcn,w);
PSM_apma=mapma(PDM,CONS,Select,typefcn,w);
PSM_hpma=mhpma(PDM,CONS,Select,typefcn,w);

if is_quality==0
    if feasibilitycheck(PSM_hpma,PDM,CONS)==false
        if feasibilitycheck(PSM_apma,PDM,CONS)==true
            PSM_hpma=PSM_apma;
        else
            if feasibilitycheck(PSM_tpma,PDM,CONS)==true
                PSM_hpma=PSM_tpma;
            end
        end
    end
else
    if feasibilitycheckq(PSM_hpma,PDM,CONS)==false
        if feasibilitycheckq(PSM_apma,PDM,CONS)==true
            PSM_hpma=PSM_apma;
        else
            if feasibilitycheckq(PSM_tpma,PDM,CONS)==true
                PSM_hpma=PSM_tpma;
            end
        end
    end
end

%% Compare results I - comparing logic networks

figure('Name','Logic networks','units','normalized','outerposition',[0 0 1 1])
subplot(1,3,1)
plot(digraph(PSM_tpma(:,1:N)));
title('Logic network of TPMa');


subplot(1,3,2)
plot(digraph(PSM_apma(:,1:N)));
title('Logic network of APMa');

subplot(1,3,3)
plot(digraph(PSM_hpma(:,1:N)));
title('Logic network of HPMa');


%% Compare results II - comparing the resource allocation

figure('Name','Resource graph - TPMa','units','normalized','outerposition',[0 0 1 1])
[BP,RESFUNC]=resfunc(PSM_tpma(:,1:N),PSM_tpma(:,end),PSM_tpma(:,N+1),PSM_tpma(:,end-nR:end-1));
drawres(BP,RESFUNC);

figure('Name','Resource graph - APMa','units','normalized','outerposition',[0 0 1 1])
[BP,RESFUNC]=resfunc(PSM_apma(:,1:N),PSM_apma(:,end),PSM_apma(:,N+1),PSM_apma(:,end-nR:end-1));
drawres(BP,RESFUNC);

figure('Name','Resource graph - HPMa','units','normalized','outerposition',[0 0 1 1])
[BP,RESFUNC]=resfunc(PSM_hpma(:,1:N),PSM_hpma(:,end),PSM_hpma(:,N+1),PSM_hpma(:,end-nR:end-1));
drawres(BP,RESFUNC);

%% Compare results III - Scheduling performances

TPT_tpma=tptsst(PSM_tpma(:,1:N),PSM_tpma(:,N+1),PSM_tpma(:,end));
TPT_apma=tptsst(PSM_apma(:,1:N),PSM_apma(:,N+1),PSM_apma(:,end));
TPT_hpma=tptsst(PSM_hpma(:,1:N),PSM_hpma(:,N+1),PSM_hpma(:,end));

TPC_tpma=tpcfast(PSM_tpma(:,1:N),PSM_tpma(:,N+2));
TPC_apma=tpcfast(PSM_apma(:,1:N),PSM_apma(:,N+2));
TPC_hpma=tpcfast(PSM_hpma(:,1:N),PSM_hpma(:,N+2));

TPR_tpma=maxresfun(PSM_tpma(:,end),PSM_tpma(:,1:N),PSM_tpma(:,N+1),PSM_tpma(:,end-nR:end-1))';
TPR_apma=maxresfun(PSM_apma(:,end),PSM_apma(:,1:N),PSM_apma(:,N+1),PSM_apma(:,end-nR:end-1))';
TPR_hpma=maxresfun(PSM_hpma(:,end),PSM_hpma(:,1:N),PSM_hpma(:,N+1),PSM_hpma(:,end-nR:end-1))';

if is_quality==1
   TPQ_tpma=tpqfast(PSM_tpma(:,1:N),PDM(:,1:N),PSM_tpma(:,N+3)); 
   TPQ_apma=tpqfast(PSM_apma(:,1:N),PDM(:,1:N),PSM_apma(:,N+3)); 
   TPQ_hpma=tpqfast(PSM_hpma(:,1:N),PDM(:,1:N),PSM_hpma(:,N+3));    
end
   
TPS_tpma=maxscore_PEM(PSM_tpma(:,1:N),PDM(:,1:N),1-PDM(:,1:N));
TPS_apma=maxscore_PEM(PSM_apma(:,1:N),PDM(:,1:N),1-PDM(:,1:N));
TPS_hpma=maxscore_PEM(PSM_hpma(:,1:N),PDM(:,1:N),1-PDM(:,1:N));

if is_quality==0
    TPMa=[1-(TPT_tpma-TPTmin)/(Ct-TPTmin),1-(TPC_tpma-TPCmin)/(Cc-TPCmin),TPR_tpma./Cr,(TPS_tpma-Cs)/(TPSmax-Cs)];
    APMa=[1-(TPT_apma-TPTmin)/(Ct-TPTmin),1-(TPC_apma-TPCmin)/(Cc-TPCmin),TPR_apma./Cr,(TPS_apma-Cs)/(TPSmax-Cs)];
    HPMa=[1-(TPT_hpma-TPTmin)/(Ct-TPTmin),1-(TPC_hpma-TPCmin)/(Cc-TPCmin),TPR_hpma./Cr,(TPS_hpma-Cs)/(TPSmax-Cs)];
else
    TPMa=[1-(TPT_tpma-TPTmin)/(Ct-TPTmin),1-(TPC_tpma-TPCmin)/(Cc-TPCmin),(TPQ_tpma-Cq)/(TPQmax-Cq),TPR_tpma./Cr,(TPS_tpma-Cs)/(TPSmax-Cs)];
    APMa=[1-(TPT_apma-TPTmin)/(Ct-TPTmin),1-(TPC_apma-TPCmin)/(Cc-TPCmin),(TPQ_apma-Cq)/(TPQmax-Cq),TPR_apma./Cr,(TPS_apma-Cs)/(TPSmax-Cs)];
    HPMa=[1-(TPT_hpma-TPTmin)/(Ct-TPTmin),1-(TPC_hpma-TPCmin)/(Cc-TPCmin),(TPQ_hpma-Cq)/(TPQmax-Cq),TPR_hpma./Cr,(TPS_hpma-Cs)/(TPSmax-Cs)];
end

xPMa=[TPMa;APMa;HPMa]'*100;

c=categorical({});
c(1)='TPT%';
c(2)='TPC%';
if is_quality==0
    for i=1:nR
        c(i+2)=strcat('TPR_',num2str(i),'%');
    end
    c(nR+3)='TPS%';
else
    c(3)='TPQ%';
    for i=1:nR
        c(i+3)=strcat('TPR_',num2str(i),'%');
    end
    c(nR+4)='TPS%';
end

figure('Name','Scheduling performance','units','normalized','outerposition',[0 0 1 1])
bar(c,xPMa)
legend({'TPMa','APMa','HPMa'},'Location','northwest','NumColumns',3)
ylabel('Performance - TPX%')
xlabel('*Note: The project plan is feasible if and only if when every performance values are positive.')

%% Risk analysis

% Phase 1 : The effect of changing demands
% Phase 2 : The effect of shocks
% Phase 3 : The effect of changing structure and priorities

dist_type=0 % 0 = beta, 1 = uniform

a=-0.1 % Optimistic scenario  (at least, most likely * 0,9)
b=+0.3 % Pessimistic scenario (at most, most likely * 1,3)

if dist_type==0
    PDM1=phase1beta(PDM,a,b);
else
    PDM1=phase1(PDM,a,b);
end

p=0.10 % Probability of shocks
s=2.00 % Shock effect ratio

PDM2=phase2(PDM1,p,s);

P=0.01 % Probability of structural change
S=0.50 % Increase/decrease ratio of priorities 

%(S>0 => increase of priorities S<0 decrease of priorities

PDM3=phase3(PDM2,P,S);

%% Draw logic networks

figure('Name','Risk analysis - structural changes','units','normalized','outerposition',[0 0 1 1])
subplot(1,4,1)
LD=triu(PDM(:,1:N),1);
ld=LD;
ld(ld<1)=0;
G=digraph(LD);
g=digraph(ld);
h=plot(G);
Tasks=1:N;
highlight(h,g,'EdgeColor','r','LineWidth',1.5);
highlight(h,floor(Tasks(diag(PDM(:,1:N))==1)),'NodeColor','r')
title('Labelled graph');
subplot(1,4,1)
LD=triu(PDM(:,1:N),1);
ld=LD;
ld(ld<1)=0;
G=digraph(LD);
g=digraph(ld);
h=plot(G);
Tasks=1:N;
highlight(h,g,'EdgeColor','r','LineWidth',1.5);
highlight(h,floor(Tasks(diag(PDM(:,1:N))==1)),'NodeColor','r')
title('Original graph');

subplot(1,4,2)
LD=triu(PDM1(:,1:N),1);
ld=LD;
ld(ld<1)=0;
G=digraph(LD);
g=digraph(ld);
h=plot(G);
Tasks=1:N;
highlight(h,g,'EdgeColor','r','LineWidth',1.5);
highlight(h,floor(Tasks(diag(PDM(:,1:N))==1)),'NodeColor','r')
title('Logic network of phase 1');

subplot(1,4,3)
LD=triu(PDM2(:,1:N),1);
ld=LD;
ld(ld<1)=0;
G=digraph(LD);
g=digraph(ld);
h=plot(G);
Tasks=1:N;
highlight(h,g,'EdgeColor','r','LineWidth',1.5);
highlight(h,floor(Tasks(diag(PDM(:,1:N))==1)),'NodeColor','r')
title('Logic network of phase 2');

subplot(1,4,4)
LD=triu(PDM3(:,1:N),1);
ld=LD;
ld(ld<1)=0;
G=digraph(LD);
g=digraph(ld);
h=plot(G);
Tasks=1:N;
highlight(h,g,'EdgeColor','r','LineWidth',1.5);
highlight(h,floor(Tasks(diag(PDM(:,1:N))==1)),'NodeColor','r')
title('Logic network of phase 3');
%% Draw demands

TPTmin1=percentt(PDM1,w,0);
TPTmax1=percentt(PDM1,w,1);
TPCmin1=percentc(PDM1,w,0);
TPCmax1=percentc(PDM1,w,1);
TPSmin1=percents(PDM1,0);
TPSmax1=percents(PDM1,1);

TPTmin2=percentt(PDM2,w,0);
TPTmax2=percentt(PDM2,w,1);
TPCmin2=percentc(PDM2,w,0);
TPCmax2=percentc(PDM2,w,1);
TPSmin2=percents(PDM2,0);
TPSmax2=percents(PDM2,1);

TPTmin3=percentt(PDM3,w,0);
TPTmax3=percentt(PDM3,w,1);
TPCmin3=percentc(PDM3,w,0);
TPCmax3=percentc(PDM3,w,1);
TPSmin3=percents(PDM3,0);
TPSmax3=percents(PDM3,1);

figure('Name','Descriptive statistics of risk effects','units','normalized','outerposition',[0 0 1 1])

if is_quality==0
    subplot(1,4,1)
    bar(categorical({'Original','Phase 1','Phase 2','Phase 3'}),[[TPTmin,TPTmax];[TPTmin1,TPTmax1];[TPTmin2,TPTmax2];[TPTmin3,TPTmax3]])
    title('Duration')
    legend({'TPT_{min}','TPT_{max}'},'Location','northwest')
    
    subplot(1,4,2)
    bar(categorical({'Original','Phase 1','Phase 2','Phase 3'}),[[TPCmin,TPCmax];[TPCmin1,TPCmax1];[TPCmin2,TPCmax2];[TPCmin3,TPCmax3]])
    title('Cost demands')
    legend({'TPC_{min}','TPC_{max}'},'Location','northwest')

    subplot(1,4,3)
    
    TPRmin1=percentr(PDM1,w,0);
    TPRmax1=percentr(PDM1,w,1);
    TPR1=[TPRmin1,TPRmax1];
    
    TPRmin2=percentr(PDM2,w,0);
    TPRmax2=percentr(PDM2,w,1);
    TPR2=[TPRmin2,TPRmax2];
    
    TPRmin3=percentr(PDM3,w,0);
    TPRmax3=percentr(PDM3,w,1);
    TPR3=[TPRmin3,TPRmax3];
    c=categorical({});
    for j=1:nR 
         c(j)=strcat('Original - TPR_',num2str(j));
         c(j+nR)=strcat('Phase 1 - TPR_',num2str(j));
         c(j+2*nR)=strcat('Phase 2 - TPR_',num2str(j));         
         c(j+3*nR)=strcat('Phase 3 - TPR_',num2str(j));         
                  
    end     
    barh([TPR;TPR1;TPR2;TPR3]);
    legend({'TPR_{min}','TPR_{max}'},'Location','northeast')
    yticklabels(c)    
    title('Resource demands')
    set(gca,'YDir','reverse')
    subplot(1,4,4)
    bar(categorical({'Original','Phase 1','Phase 2','Phase 3'}),[[TPSmin,TPSmax];[TPSmin1,TPSmax1];[TPSmin2,TPSmax2];[TPSmin3,TPSmax3]])
    title('Score/Scope demands')
    legend({'TPS_{min}','TPS_{max}'},'Location','northwest')
else
    
    subplot(1,5,1)
    bar(categorical({'Original','Phase 1','Phase 2','Phase 3'}),[[TPTmin,TPTmax];[TPTmin1,TPTmax1];[TPTmin2,TPTmax2];[TPTmin3,TPTmax3]])
    title('Duration')
    legend({'TPT_{min}','TPT_{max}'},'Location','northwest')
    
    subplot(1,5,2)
    bar(categorical({'Original','Phase 1','Phase 2','Phase 3'}),[[TPCmin,TPCmax];[TPCmin1,TPCmax1];[TPCmin2,TPCmax2];[TPCmin3,TPCmax3]])
    title('Cost demands')
    legend({'TPC_{min}','TPC_{max}'},'Location','northwest')

    subplot(1,5,3)
    
    TPQmin1=percentq(PDM1,w,0);
    TPQmax1=percentq(PDM1,w,1);
    
    TPQmin2=percentq(PDM2,w,0);
    TPQmax2=percentq(PDM2,w,1);
    
    TPQmin3=percentq(PDM3,w,0);
    TPQmax3=percentq(PDM3,w,1);
    
    
    bar(categorical({'Original','Phase 1','Phase 2','Phase 3'}),[[TPQmin,TPQmax];[TPQmin1,TPQmax1];[TPQmin2,TPQmax2];[TPQmin3,TPQmax3]])
    title('Quality parameters')
    legend({'TPQ_{min}','TPQ_{max}'},'Location','northwest')

    subplot(1,5,4)
    
    TPRmin1=percentqr(PDM1,w,0);
    TPRmax1=percentqr(PDM1,w,1);
    TPR1=[TPRmin1,TPRmax1];
    
    TPRmin2=percentqr(PDM2,w,0);
    TPRmax2=percentqr(PDM2,w,1);
    TPR2=[TPRmin2,TPRmax2];
    
    TPRmin3=percentqr(PDM3,w,0);
    TPRmax3=percentqr(PDM3,w,1);
    TPR3=[TPRmin3,TPRmax3];
    c=categorical({});
    for j=1:nR 
         c(j)=strcat('Original - TPR_',num2str(j));
         c(j+nR)=strcat('Phase 1 - TPR_',num2str(j));
         c(j+2*nR)=strcat('Phase 2 - TPR_',num2str(j));         
         c(j+3*nR)=strcat('Phase 3 - TPR_',num2str(j));         
                  
    end     
    barh([TPR;TPR1;TPR2;TPR3]);
    legend({'TPR_{min}','TPR_{max}'},'Location','northeast')
    yticklabels(c)    
    title('Resource demands')
    set(gca,'YDir','reverse')
    subplot(1,5,5)
    bar(categorical({'Original','Phase 1','Phase 2','Phase 3'}),[[TPSmin,TPSmax];[TPSmin1,TPSmax1];[TPSmin2,TPSmax2];[TPSmin3,TPSmax3]])
    title('Score/Scope demands')
    legend({'TPS_{min}','TPS_{max}'},'Location','northwest')
end

%% Optimization or phases of risk analysis

PSM_tpma1=mtpma(PDM1,CONS,Select,typefcn,w);
PSM_apma1=mapma(PDM1,CONS,Select,typefcn,w);
PSM_hpma1=mhpma(PDM1,CONS,Select,typefcn,w);

if is_quality==0
    if feasibilitycheck(PSM_hpma1,PDM1,CONS)==false
        if feasibilitycheck(PSM_apma1,PDM1,CONS)==true
            PSM_hpma1=PSM_apma1;
        else
            if feasibilitycheck(PSM_tpma1,PDM1,CONS)==true
                PSM_hpma1=PSM_tpma1;
            end
        end
    end
else
    if feasibilitycheckq(PSM_hpma1,PDM1,CONS)==false
        if feasibilitycheckq(PSM_apma1,PDM1,CONS)==true
            PSM_hpma1=PSM_apma1;
        else
            if feasibilitycheckq(PSM_tpma1,PDM1,CONS)==true
                PSM_hpma1=PSM_tpma1;
            end
        end
    end
end

PSM_tpma2=mtpma(PDM2,CONS,Select,typefcn,w);
PSM_apma2=mapma(PDM2,CONS,Select,typefcn,w);
PSM_hpma2=mhpma(PDM2,CONS,Select,typefcn,w);

if is_quality==0
    if feasibilitycheck(PSM_hpma2,PDM2,CONS)==false
        if feasibilitycheck(PSM_apma2,PDM2,CONS)==true
            PSM_hpma2=PSM_apma2;
        else
            if feasibilitycheck(PSM_tpma2,PDM2,CONS)==true
                PSM_hpma2=PSM_tpma2;
            end
        end
    end
else
    if feasibilitycheckq(PSM_hpma2,PDM2,CONS)==false
        if feasibilitycheckq(PSM_apma2,PDM2,CONS)==true
            PSM_hpma2=PSM_apma2;
        else
            if feasibilitycheckq(PSM_tpma2,PDM2,CONS)==true
                PSM_hpma2=PSM_tpma2;
            end
        end
    end
end

PSM_tpma3=mtpma(PDM3,CONS,Select,typefcn,w);
PSM_apma3=mapma(PDM3,CONS,Select,typefcn,w);
PSM_hpma3=mhpma(PDM3,CONS,Select,typefcn,w);

if is_quality==0
    if feasibilitycheck(PSM_hpma3,PDM3,CONS)==false
        if feasibilitycheck(PSM_apma3,PDM3,CONS)==true
            PSM_hpma3=PSM_apma3;
        else
            if feasibilitycheck(PSM_tpma3,PDM3,CONS)==true
                PSM_hpma3=PSM_tpma3;
            end
        end
    end
else
    if feasibilitycheckq(PSM_hpma3,PDM3,CONS)==false
        if feasibilitycheckq(PSM_apma3,PDM3,CONS)==true
            PSM_hpma3=PSM_apma3;
        else
            if feasibilitycheckq(PSM_tpma3,PDM3,CONS)==true
                PSM_hpma3=PSM_tpma3;
            end
        end
    end
end

%% Compare results  of risk analysis I - comparing logic networks

figure('Name','Logic networks of results of risk analysis','units','normalized','outerposition',[0 0 1 1])
subplot(4,3,1)
plot(digraph(PSM_tpma(:,1:N)));
title('Original network of TPMa');

subplot(4,3,2)
plot(digraph(PSM_apma(:,1:N)));
title('Original network of APMa');

subplot(4,3,3)
plot(digraph(PSM_hpma(:,1:N)));
title('Original network of HPMa');

subplot(4,3,4)
plot(digraph(PSM_tpma1(:,1:N)));
title('Phase 1 network of TPMa');

subplot(4,3,5)
plot(digraph(PSM_apma1(:,1:N)));
title('Phase 1 network of APMa');

subplot(4,3,6)
plot(digraph(PSM_hpma1(:,1:N)));
title('Phase 1 network of HPMa');

subplot(4,3,7)
plot(digraph(PSM_tpma2(:,1:N)));
title('Phase 2 network of TPMa');

subplot(4,3,8)
plot(digraph(PSM_apma2(:,1:N)));
title('Phase 2 network of APMa');

subplot(4,3,9)
plot(digraph(PSM_hpma2(:,1:N)));
title('Phase 2 network of HPMa');

subplot(4,3,10)
plot(digraph(PSM_tpma3(:,1:N)));
title('Phase 3 network of TPMa');

subplot(4,3,11)
plot(digraph(PSM_apma3(:,1:N)));
title('Phase 3 network of APMa');

subplot(4,3,12)
plot(digraph(PSM_hpma3(:,1:N)));
title('Phase 3 network of HPMa');

%% Compare results of risk analysis II - Scheduling performances

TPT_tpma1=tptsst(PSM_tpma1(:,1:N),PSM_tpma1(:,N+1),PSM_tpma1(:,end));
TPT_apma1=tptsst(PSM_apma1(:,1:N),PSM_apma1(:,N+1),PSM_apma1(:,end));
TPT_hpma1=tptsst(PSM_hpma1(:,1:N),PSM_hpma1(:,N+1),PSM_hpma1(:,end));

TPC_tpma1=tpcfast(PSM_tpma1(:,1:N),PSM_tpma1(:,N+2));
TPC_apma1=tpcfast(PSM_apma1(:,1:N),PSM_apma1(:,N+2));
TPC_hpma1=tpcfast(PSM_hpma1(:,1:N),PSM_hpma1(:,N+2));

TPR_tpma1=maxresfun(PSM_tpma1(:,end),PSM_tpma1(:,1:N),PSM_tpma1(:,N+1),PSM_tpma1(:,end-nR:end-1))';
TPR_apma1=maxresfun(PSM_apma1(:,end),PSM_apma1(:,1:N),PSM_apma1(:,N+1),PSM_apma1(:,end-nR:end-1))';
TPR_hpma1=maxresfun(PSM_hpma1(:,end),PSM_hpma1(:,1:N),PSM_hpma1(:,N+1),PSM_hpma1(:,end-nR:end-1))';

if is_quality==1
   TPQ_tpma1=tpqfast(PSM_tpma1(:,1:N),PDM(:,1:N),PSM_tpma1(:,N+3)); 
   TPQ_apma1=tpqfast(PSM_apma1(:,1:N),PDM(:,1:N),PSM_apma1(:,N+3)); 
   TPQ_hpma1=tpqfast(PSM_hpma1(:,1:N),PDM(:,1:N),PSM_hpma1(:,N+3));    
end
   
TPS_tpma1=maxscore_PEM(PSM_tpma1(:,1:N),PDM(:,1:N),1-PDM(:,1:N));
TPS_apma1=maxscore_PEM(PSM_apma1(:,1:N),PDM(:,1:N),1-PDM(:,1:N));
TPS_hpma1=maxscore_PEM(PSM_hpma1(:,1:N),PDM(:,1:N),1-PDM(:,1:N));

if is_quality==0
    TPMa1=[1-(TPT_tpma1-TPTmin1)/(Ct-TPTmin1),1-(TPC_tpma1-TPCmin1)/(Cc-TPCmin1),TPR_tpma1./Cr,(TPS_tpma1-Cs)/(TPSmax1-Cs)];
    APMa1=[1-(TPT_apma1-TPTmin1)/(Ct-TPTmin1),1-(TPC_apma1-TPCmin1)/(Cc-TPCmin1),TPR_apma1./Cr,(TPS_apma1-Cs)/(TPSmax1-Cs)];
    HPMa1=[1-(TPT_hpma1-TPTmin1)/(Ct-TPTmin1),1-(TPC_hpma1-TPCmin1)/(Cc-TPCmin1),TPR_hpma1./Cr,(TPS_hpma1-Cs)/(TPSmax1-Cs)];
else
    TPMa1=[1-(TPT_tpma1-TPTmin1)/(Ct-TPTmin1),1-(TPC_tpma1-TPCmin1)/(Cc-TPCmin1),(TPQ_tpma1-Cq)/(TPQmax1-Cq),TPR_tpma1./Cr,(TPS_tpma1-Cs)/(TPSmax1-Cs)];
    APMa1=[1-(TPT_apma1-TPTmin1)/(Ct-TPTmin1),1-(TPC_apma1-TPCmin1)/(Cc-TPCmin1),(TPQ_apma1-Cq)/(TPQmax1-Cq),TPR_apma1./Cr,(TPS_apma1-Cs)/(TPSmax1-Cs)];
    HPMa1=[1-(TPT_hpma1-TPTmin1)/(Ct-TPTmin1),1-(TPC_hpma1-TPCmin1)/(Cc-TPCmin1),(TPQ_hpma1-Cq)/(TPQmax1-Cq),TPR_hpma1./Cr,(TPS_hpma1-Cs)/(TPSmax1-Cs)];
end

xPMa1=[TPMa1;APMa1;HPMa1]'*100;

TPT_tpma2=tptsst(PSM_tpma2(:,1:N),PSM_tpma2(:,N+1),PSM_tpma2(:,end));
TPT_apma2=tptsst(PSM_apma2(:,1:N),PSM_apma2(:,N+1),PSM_apma2(:,end));
TPT_hpma2=tptsst(PSM_hpma2(:,1:N),PSM_hpma2(:,N+1),PSM_hpma2(:,end));

TPC_tpma2=tpcfast(PSM_tpma2(:,1:N),PSM_tpma2(:,N+2));
TPC_apma2=tpcfast(PSM_apma2(:,1:N),PSM_apma2(:,N+2));
TPC_hpma2=tpcfast(PSM_hpma2(:,1:N),PSM_hpma2(:,N+2));

TPR_tpma2=maxresfun(PSM_tpma2(:,end),PSM_tpma2(:,1:N),PSM_tpma2(:,N+1),PSM_tpma2(:,end-nR:end-1))';
TPR_apma2=maxresfun(PSM_apma2(:,end),PSM_apma2(:,1:N),PSM_apma2(:,N+1),PSM_apma2(:,end-nR:end-1))';
TPR_hpma2=maxresfun(PSM_hpma2(:,end),PSM_hpma2(:,1:N),PSM_hpma2(:,N+1),PSM_hpma2(:,end-nR:end-1))';

if is_quality==1
   TPQ_tpma2=tpqfast(PSM_tpma2(:,1:N),PDM(:,1:N),PSM_tpma2(:,N+3)); 
   TPQ_apma2=tpqfast(PSM_apma2(:,1:N),PDM(:,1:N),PSM_apma2(:,N+3)); 
   TPQ_hpma2=tpqfast(PSM_hpma2(:,1:N),PDM(:,1:N),PSM_hpma2(:,N+3));    
end
   
TPS_tpma2=maxscore_PEM(PSM_tpma2(:,1:N),PDM(:,1:N),1-PDM(:,1:N));
TPS_apma2=maxscore_PEM(PSM_apma2(:,1:N),PDM(:,1:N),1-PDM(:,1:N));
TPS_hpma2=maxscore_PEM(PSM_hpma2(:,1:N),PDM(:,1:N),1-PDM(:,1:N));

if is_quality==0
    TPMa2=[1-(TPT_tpma2-TPTmin2)/(Ct-TPTmin2),1-(TPC_tpma2-TPCmin2)/(Cc-TPCmin2),TPR_tpma2./Cr,(TPS_tpma2-Cs)/(TPSmax2-Cs)];
    APMa2=[1-(TPT_apma2-TPTmin2)/(Ct-TPTmin2),1-(TPC_apma2-TPCmin2)/(Cc-TPCmin2),TPR_apma2./Cr,(TPS_apma2-Cs)/(TPSmax2-Cs)];
    HPMa2=[1-(TPT_hpma2-TPTmin2)/(Ct-TPTmin2),1-(TPC_hpma2-TPCmin2)/(Cc-TPCmin2),TPR_hpma2./Cr,(TPS_hpma2-Cs)/(TPSmax2-Cs)];
else
    TPMa2=[1-(TPT_tpma2-TPTmin2)/(Ct-TPTmin2),1-(TPC_tpma2-TPCmin2)/(Cc-TPCmin2),(TPQ_tpma2-Cq)/(TPQmax2-Cq),TPR_tpma2./Cr,(TPS_tpma2-Cs)/(TPSmax2-Cs)];
    APMa2=[1-(TPT_apma2-TPTmin2)/(Ct-TPTmin2),1-(TPC_apma2-TPCmin2)/(Cc-TPCmin2),(TPQ_apma2-Cq)/(TPQmax2-Cq),TPR_apma2./Cr,(TPS_apma2-Cs)/(TPSmax2-Cs)];
    HPMa2=[1-(TPT_hpma2-TPTmin2)/(Ct-TPTmin2),1-(TPC_hpma2-TPCmin2)/(Cc-TPCmin2),(TPQ_hpma2-Cq)/(TPQmax2-Cq),TPR_hpma2./Cr,(TPS_hpma2-Cs)/(TPSmax2-Cs)];
end

xPMa2=[TPMa2;APMa2;HPMa2]'*100;

TPT_tpma3=tptsst(PSM_tpma3(:,1:N),PSM_tpma3(:,N+1),PSM_tpma3(:,end));
TPT_apma3=tptsst(PSM_apma3(:,1:N),PSM_apma3(:,N+1),PSM_apma3(:,end));
TPT_hpma3=tptsst(PSM_hpma3(:,1:N),PSM_hpma3(:,N+1),PSM_hpma3(:,end));

TPC_tpma3=tpcfast(PSM_tpma3(:,1:N),PSM_tpma3(:,N+2));
TPC_apma3=tpcfast(PSM_apma3(:,1:N),PSM_apma3(:,N+2));
TPC_hpma3=tpcfast(PSM_hpma3(:,1:N),PSM_hpma3(:,N+2));

TPR_tpma3=maxresfun(PSM_tpma3(:,end),PSM_tpma3(:,1:N),PSM_tpma3(:,N+1),PSM_tpma3(:,end-nR:end-1))';
TPR_apma3=maxresfun(PSM_apma3(:,end),PSM_apma3(:,1:N),PSM_apma3(:,N+1),PSM_apma3(:,end-nR:end-1))';
TPR_hpma3=maxresfun(PSM_hpma3(:,end),PSM_hpma3(:,1:N),PSM_hpma3(:,N+1),PSM_hpma3(:,end-nR:end-1))';

if is_quality==1
   TPQ_tpma3=tpqfast(PSM_tpma3(:,1:N),PDM(:,1:N),PSM_tpma3(:,N+3)); 
   TPQ_apma3=tpqfast(PSM_apma3(:,1:N),PDM(:,1:N),PSM_apma3(:,N+3)); 
   TPQ_hpma3=tpqfast(PSM_hpma3(:,1:N),PDM(:,1:N),PSM_hpma3(:,N+3));    
end
   
TPS_tpma3=maxscore_PEM(PSM_tpma3(:,1:N),PDM(:,1:N),1-PDM(:,1:N));
TPS_apma3=maxscore_PEM(PSM_apma3(:,1:N),PDM(:,1:N),1-PDM(:,1:N));
TPS_hpma3=maxscore_PEM(PSM_hpma3(:,1:N),PDM(:,1:N),1-PDM(:,1:N));

if is_quality==0
    TPMa3=[1-(TPT_tpma3-TPTmin3)/(Ct-TPTmin3),1-(TPC_tpma3-TPCmin3)/(Cc-TPCmin3),TPR_tpma3./Cr,(TPS_tpma3-Cs)/(TPSmax3-Cs)];
    APMa3=[1-(TPT_apma3-TPTmin3)/(Ct-TPTmin3),1-(TPC_apma3-TPCmin3)/(Cc-TPCmin3),TPR_apma3./Cr,(TPS_apma3-Cs)/(TPSmax3-Cs)];
    HPMa3=[1-(TPT_hpma3-TPTmin3)/(Ct-TPTmin3),1-(TPC_hpma3-TPCmin3)/(Cc-TPCmin3),TPR_hpma3./Cr,(TPS_hpma3-Cs)/(TPSmax3-Cs)];
else
    TPMa3=[1-(TPT_tpma3-TPTmin3)/(Ct-TPTmin3),1-(TPC_tpma3-TPCmin3)/(Cc-TPCmin3),(TPQ_tpma3-Cq)/(TPQmax3-Cq),TPR_tpma3./Cr,(TPS_tpma3-Cs)/(TPSmax3-Cs)];
    APMa3=[1-(TPT_apma3-TPTmin3)/(Ct-TPTmin3),1-(TPC_apma3-TPCmin3)/(Cc-TPCmin3),(TPQ_apma3-Cq)/(TPQmax3-Cq),TPR_apma3./Cr,(TPS_apma3-Cs)/(TPSmax3-Cs)];
    HPMa3=[1-(TPT_hpma3-TPTmin3)/(Ct-TPTmin3),1-(TPC_hpma3-TPCmin3)/(Cc-TPCmin3),(TPQ_hpma3-Cq)/(TPQmax3-Cq),TPR_hpma3./Cr,(TPS_hpma3-Cs)/(TPSmax3-Cs)];
end

xPMa3=[TPMa3;APMa3;HPMa3]'*100;

figure('Name','Scheduling performance after risk analysis','units','normalized','outerposition',[0 0 1 1])
bar3([xPMa,xPMa1,xPMa2,xPMa3])
view(-15,15)
c=categorical({});
if is_quality==0
    c(1)='TPT%';
    c(2)='TPC%';
    for i=1:nR
        c(2+i)=strcat('TPR_',num2str(i),'%');
    end
    c(3+nR)='TPS%';
else
    c(1)='TPT%';
    c(2)='TPC%';
    c(3)='TPQ%';
    for i=1:nR
        c(3+i)=strcat('TPR_',num2str(i),'%');
    end
    c(4+nR)='TPS%';
end
yticklabels(c);
xticklabels({'Original - TPMa','Original - APMa','Original - HPMa','Phase 1 - TPMa','Phase 1 - APMa','Phase 1 - HPMa','Phase 2 - TPMa','Phase 2 - APMa','Phase 2 - HPMa','Phase 3 - TPMa','Phase 3 - APMa','Phase 3 - HPMa'});
title('Scheduling performances after risk analysis')
