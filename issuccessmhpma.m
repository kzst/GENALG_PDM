%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Feasibility check of hybrid project management agent (APMa)
%----------------
%Output:
%Flag: %0=failed (infeasible for both sCONS and cCONS); 1=challanged 
%(feasible for cCONS but infeasible for sCONS); 2=success (feasible for
%sCONS)
%---------------- 
%Inputs:
%PDM: N by N+(M-N)*W matrix of the stochastic project plan. 
%PDM=[PEM,TD,CD{,QD}{,RD}], where PEM is an N by N upper triangular matrix 
%of logic domain, 
% TD is an N by w matrix of task durations,
% CD is an N by w matrix of cost demands,
%depends on the problem selection
% QD is an N by w matrix of quality parameters,
%depends on the resource allocation is performed or not
% RD is an N by w*nR matrix of resource demands.
%sCONS,cCONS: 3..3+nR+1 (depends on the problem selection) elements row vector of 
%[Ct,Cc,{Cq,CR,}Cs] values, where
% Ct is the time constraint (max constraint, required in all problems)
% Cc is the cost constraint (max constraint, required in all problems)
% Cq is the quality constraint (min constraint, required in TQCTP problems))
% CR are the resource constraints (max constraint, required in 
%resource-allocation problems)
% Cs is the score constraint (min constraint, required in all problems)
%Select: 1-18 scalar of problem selection
% 1: APSP: Agile Project Scheduling Problem
% 2: APSPq: Agile Project Scheduling Problem with quality parameters
% 3: HDTCTP: Hybrid Discrete Time-Cost Trade-off Problem
% 4: HCTCTP: Hybrid Continouos Time-Cost Trade-off Problem
% 5: HDTQCTP: Hybrid Discrete Time-Quality-Cost Trade-off Problem
% 6: HCTQCTP: Hybrid Continouos Time-Quality-Cost Trade-off Problem
% 7: RC-APSP: Resource-Constraint Agile Project Scheduling Problem 
% 8: PO-APSP: Pareto-Optimal resource-allocated Agile Project Scheduling 
%Problem 
% 9: RC-HDTCTP: Resource-Constraint Hybrid Discrete Time-Cost Trade-off 
%Problem
%10: RC-HCTCTP: Resource-Constraint Hybrid Continouos Time-Cost Trade-off 
%Problem
%11: PO-HDTCTP: Pareto-Optimal resource-allocated Hybrid Discrete Time-Cost 
%Trade-off Problem
%12: PO-HCTCTP: Pareto-Optimal resource-allocated Hybrid Continouos 
%Time-Cost Trade-off Problem
%13: RC-APSPq: Resource-Constraint Agile Project Scheduling Problem 
%with quality parameters
%14: PO-APSPq: Pareto-Optimal resource-allocated Agile Project Scheduling 
%with quality parameters
%15: RC-HDTCTP: Resource-Constraint Hybrid Discrete Time-Quality-Cost 
%Trade-off Problem
%16: PO-HDTCTP: Pareto-Optimal resource-allocated Hybrid Discrete 
%Time-Quality-Cost Trade-off Problem
%17: RC-HCTCTP: Resource-Constraint Hybrid Continouos Time-Quality-Cost 
%Trade-off Problem
%18: PO-HCTCTP: Pareto-Optimal resource-allocated Hybrid Continouos 
%Time-Quality-Cost Trade-off Problem
%typetfcn: Type of target function 
% {0=maxTPQ,} 1=minTPT, 2=minTPC, 3=maxTPS, {4=minUF,} ~ composite
%w: Number of modes
%---------------- 
%Usage:
%Flag=issuccessmhpma(PDM,sCONS,cCONS,Select,typefcn,w)
%---------------- 
%Example:
%PDM=generatepdmq(30,0.05,0,20,30,20,2,2,2);
%Select=13;
%typefcn=999; %let target function be a composite target function
%w=3;
%Flag=issuccessmhpma(PDM,sCONS,cCONS,Select,typefcn,w)
%---------------- 
%Prepositions and Requirements:
%1.)The logic domian of the input PDM matrix must be an upper triangular 
%matrix, where the matrix elements are between 0 to 1 interval.
%2.)Before running the evaluation, all matrix element (PEMi,j)of the 
%original logic domain is converted to 0, if the PEMi,j<0.5, otherwise this
%element is converted to 1. It merans, if the score of inclusion (Pi,j) is 
%greather than the score of exclusion (Qi,j). In this case the score value
%of the project scenario for converted matrix is highest.
%3.)The number of modes(w) is a positive integer. In case of countinous
%trade-off problems w=2 and in case of selections the w shuld be 1.
%4.)The number of resources (nR) is a positive number
%5.)The elements of Time/Cost/Quality/Resource Domains are positive real 
%numbers.
%6.)Usually a monotonity are assumed, which, means: if tk,i<tk,j => 
%ck,i>=ck,j and qk,i<=qk,j, where 1<=k<=N is the k-th task, 1<=i,j<=w are 
%the selected modes. This asumptation is only used in continous trade-off
%problems. 
%7.)In the case of TPMa there is no flexible task dependency and uncertain 
%task completions (see 2.), however, multiple modes are considered.

function Flag=issuccessmhpma(PDM,sCONS,cCONS,Select,typefcn,w)
Flag=0; %0=failed;1=challanged;2=success
PEM=PDM(:,1:size(PDM,1)); % N by N matrix of logic domain
P=PEM; % Initial N by N score matrix of task/dependency inclusion
Q=ones(size(PEM,1))-P; %Q=1-P
sCR=[]; %Resource constraints
if ((Select==2)||(Select==5)||(Select==6)||(Select>=13)) %If quality 
    if numel(sCONS)>4                      %paramters are considered
        sCR=sCONS(4:end-1)';               %[Ct,Cc,Cq,CR',Cs]
    end
else
    if numel(sCONS)>3                      %[Ct,Cc,CR',Cs]
        sCR=sCONS(3:end-1)';
    end
end
PSM=mhpma(PDM,sCONS,Select,typefcn,w);     %Run PMa
N=size(PSM,1);
M=size(PSM,2);
DSM=PSM(:,1:size(PSM,1));
TD=PSM(:,size(PSM,1)+1);
CD=PSM(:,size(PSM,1)+2);
QD=[];
if ((Select==2)||(Select==5)||(Select==6)||(Select>=13))
    QD=PSM(:,size(PSM,1)+3);
end
SST=PSM(:,end);
RD=[];
if ((Select==2)||(Select==5)||(Select==6)||(Select>=13))
    if M>N+4
        RD=PSM(:,size(PSM,1)+4:end-1);
    end
else
    if M>N+3
        RD=PSM(:,size(PSM,1)+3:end-1);
    end
end
scons=sCONS;           %If quality parameter is considered...
if ((Select==2)||(Select==5)||(Select==6)||(Select>=13)) 
    if numel(sCR)>0          %If resources are considered
        scons(3)=-sCONS(3);     %Cq is maximal constraint
        scons(end)=-sCONS(end); %Cs is maximal constraint
        sC=min(scons-[tptfast(DSM,TD),tpcfast(DSM,CD),...
            -tpqfast(DSM,P,QD),maxresfun(SST,DSM,TD,RD)',...
            -maxscore_PEM(DSM,P,Q)]);
    else
        scons(3)=-sCONS(3);     %Cq is maximal constraint
        scons(end)=-sCONS(end); %Cs is maximal constraint
        sC=min(scons-[tptfast(DSM,TD),tpcfast(DSM,CD),...
            -tpqfast(DSM,P,QD),...
            -maxscore_PEM(DSM,P,Q)]);
    end
else
    if numel(sCR)>0          %If resources are considered
        scons(end)=-sCONS(end); %Cs is maximal constraint
        sC=min(scons-[tptfast(DSM,TD),tpcfast(DSM,CD),...
            maxresfun(SST,DSM,TD,RD)',-maxscore_PEM(DSM,P,Q)]);
    else
        scons(end)=-sCONS(end); %Cs is maximal constraint
        sC=min(scons-[tptfast(DSM,TD),tpcfast(DSM,CD),...
            -maxscore_PEM(DSM,P,Q)]);
    end
end
if (sC>=0)
    Flag=2; %Success
else
    PSM=mhpma(PDM,cCONS,Select,typefcn,w); %Run PMa for relaxed constraint
    N=size(PSM,1);
    M=size(PSM,2);
    DSM=PSM(:,1:size(PSM,1)); 
    TD=PSM(:,size(PSM,1)+1);  
    CD=PSM(:,size(PSM,1)+2);
    QD=[];
    if ((Select==2)||(Select==5)||(Select==6)||(Select>=13))
        QD=PSM(:,size(PSM,1)+3);
    end
    SST=PSM(:,end);
    RD=[];
    if ((Select==2)||(Select==5)||(Select==6)||(Select>=13))
        if M>N+4
            RD=PSM(:,size(PSM,1)+4:end-1);
        end
    else
        if M>N+3
            RD=PSM(:,size(PSM,1)+3:end-1);
        end
    end
    ccons=cCONS;        %If quality parameter is considered...
    if ((Select==2)||(Select==5)||(Select==6)||(Select>=13))
        if numel(sCR)>0          %If resources are considered
            ccons(3)=-cCONS(3);     %Cq is maximal constraint
            ccons(end)=-cCONS(end); %Cs is maximal constraint
            cC=min(ccons-[tptfast(DSM,TD),tpcfast(DSM,CD),...
                -tpqfast(DSM,P,QD),maxresfun(SST,DSM,TD,RD)',...
                -maxscore_PEM(DSM,P,Q)]);
        else
            ccons(3)=-cCONS(3);     %Cq is maximal constraint
            ccons(end)=-cCONS(end); %Cs is maximal constraint
            cC=min(ccons-[tptfast(DSM,TD),tpcfast(DSM,CD),...
                -tpqfast(DSM,P,QD),...
                -maxscore_PEM(DSM,P,Q)]);
        end
    else
        if numel(sCR)>0
            ccons(end)=-cCONS(end); %Cs is maximal constraint
            cC=min(ccons-...       
                [tptfast(DSM,TD),tpcfast(DSM,CD),...
                maxresfun(SST,DSM,TD,RD)',...
                -maxscore_PEM(DSM,P,Q)]);
        else
            ccons(end)=-cCONS(end); %Cs is maximal constraint
            cC=min(ccons-...
                [tptfast(DSM,TD),tpcfast(DSM,CD),...
                -maxscore_PEM(DSM,P,Q)]);
        end
    end
    if (cC>=0)
        Flag=1; %Challenged
    end
end