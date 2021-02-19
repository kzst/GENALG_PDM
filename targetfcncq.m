%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%GA target function (fitness function) for
%Hybrid Continouos Time-Quality-Cost Trade-off Problems (HCTQCTP)
%----------------
%Standard output:
%target: is a real value, which has to be minimized
%---------------- 
%Standard inputs:
%chromosome: is a set of genomes
%---------------- 
%Usage:
%target=targetfcncq(chromosome)
%---------------- 
%Example:
%This is is not executed as a stand-alone fuinction. Instead of specifying
%global values this function is cannot be run.
%---------------- 
%Prepositions and Requirements:
%1.)Global variables must be specified.

function target=targetfcncq(chromosome)
global P Q T C q Ct Cc Cq Cs Typ QD
%global variables: 
% P is the score matrix of inclusion task dependencies/completions 
%the logic domain (initially PEM=P)
% Q is the score matrix of exclusion task dependencies/completions. 
%In this simulation Q=1-P
% T is the column vector of time demands (=the time domain)
% C is the column vector of cost demands (=the cost domain)
% q is the column vector of quality parameters (=the quality domain)
% Ct is the time constraint
% Cc is the cost constraint
% Cq is the quality constraint
% Cs is the score constraint
% Typ is the selection number of the target functions (see above)

%-----------Calculation of maximal, minimal values--------------%
%-----------for determining the effectiveness of GA-------------%

dsm=floor(P);
dsmdiag=diag(dsm);
if size(T,2)>1
    t=min(T')';
else
    t=T;
end
if size(C,2)>1
    c=min(C')';
else
    c=C;
end
if size(q,2)>1
    qq=max(q')';
else
    qq=q;
end
t(dsmdiag==0)=0;
c(dsmdiag==0)=0;
qq(dsmdiag==0)=0;
dT=1;
dC=1;
dS=1;
dQ=1;
TPTmin=tptfast(dsm,t);
TPCmin=tpcfast(dsm,c);
TPSmax=maxscore_PEM(round(P),P,Q);
TPQmax=tpqfastqd(round(P),P,QD,qq);
if (Ct-TPTmin)>0
    dT=Ct-TPTmin;
end
if (Cc-TPCmin)>0
    dC=Cc-TPCmin;
end
if (TPSmax-Cs)>0
    dS=TPSmax-Cs;
end
if (TPQmax-Cq)>0
    dQ=TPQmax-Cq;
end

%----End of calculating minimal, maximal values----

PopSize=numel(chromosome(:,1)); %Number of chromosomes
nA=numel(P(P>0&P<1)); %Number of uncertain activities/dependencies
N=size(P,1); %Number of activities
chromosome(:,1:nA)=round(chromosome(:,1:nA));
for i=1:PopSize
    s(i).f1 = chromosome(i,:);
end

tpts=zeros(PopSize,1); %PopSize by 1 vector of TPTs
tpcs=zeros(PopSize,1); %PopSize by 1 vector of TPCs
tpqs=zeros(PopSize,1); %PopSize by 1 vector of TPQs
c=zeros(1,PopSize); %Penalty part of target function
d=zeros(PopSize,1); %Rewarding part of target function

A=arrayfun(@(x) updatepemc(x.f1), s,'UniformOutput',false); %set of PSMs
for i=1:PopSize
    PSM=A{i};
    DSM=round(PSM(:,1:size(PSM,1)));
    MODES=PSM(:,end);
    TD=zeros(N,1);
    CD=zeros(N,1);
    qd=zeros(N,1);
    for I=1:N
        if MODES(I)==0
            TD(I)=0;
            CD(I)=0;
            qd(i)=0;
        else
            TD(I)=MODES(I);
            CD(I)=timetocost(TD(I),T(I,:),C(I,:));
            qd(I)=timetoquality(TD(I),T(I,:),q(I,:));
        end
    end
    tpts(i)=tptfast(DSM,TD);
    tpqs(i)=tpqfastqd(DSM,P,QD,qd);
    tpcs(i)=diag(DSM)'*CD;
    c(i)=max([tpts(i)-Ct,tpcs(i)-Cc,Cs-maxscore_PEM(DSM,P,Q),Cq-tpqs(i)]);
    if c(i)>0 %Penalty case: if constraints are NOT satisfied
        c(i)=1e6*(exp(c(i))-1);
        d(i)=c(i);
    else %Rewarding case: if constraints are satisfied
        d(i)=1-(((Ct-tpts(i))/dT)^(1/4)*((Cc-tpcs(i))/dC)^(1/4)*...
           ((maxscore_PEM(A{i},P,Q)-Cs)/dS)^(1/4)*((tpqs(i)-Cq)/dQ)^(1/4));
        if Typ<3
            c(i)=0; %Penalty part is deleted in the rewarding cases
        end
    end
end
typ=Typ;
switch typ
    case 0 %maximize TPQ
        target=-tpqs+c';
    case 1 %minimize TPT
        target=tpts+c';
    case 2 %minimize TPC
        target=tpcs+c';
    case 3 %minimize TPS
        target=-(maxscore_PEM_cell(A,P,Q)'*100+maxscore_SNPM_cell(A,P,Q)');
    otherwise %maximize efectiveness
       target=d;
end