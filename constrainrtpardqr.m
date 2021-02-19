%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Calculate the constraints for the resource-constraint hybrid 
%discrete time-quality-cost trade-off problems
%----------------
%Outputs:
%c is vector of inequities.
%ceq is null vector of Eq.s
%---------------- 
%Input:
%chromosome: set of populations 
%---------------- 
%Usage:
%[c,ceq]=constrainrtpardqr(chromosome)
%---------------- 
%Example:
%This is is not executed as a stand-alone fuinction. Instead of specifying
%global values this function is cannot be run.
%---------------- 
%Prepositions and Requirements:
%1.)Global values P,Q,T,C,R,q,Ct,Cc,CR,Cq,Cs must be specified.
%2.)P should be an upper triangular matrix of the logic domain.

function [c,ceq]=constrainrtpardqr(chromosome)
global P Q T C R q Ct Cc CR Cq Cs

PopSize=numel(chromosome(:,1)); %The size of the population
N=size(P,1); %Number of activities
w=size(T,2); %Number of modes
nR=size(R,2)/w; %Number of resources
c=zeros(PopSize,3+nR); %Initialize the PopSize by 3+nR matrix of inequities
tpts=zeros(PopSize,1); %Initialize the vector of TPTs
tpcs=zeros(PopSize,1); %Initialize the vector of TPCs
tpqs=zeros(PopSize,1); %Initialize the vector of TPQs

for I=1:PopSize
    PSM=updatepemdr(chromosome(I,:));%PSM=[DSM,MODES,SST] = Logic Domain, 
                                 %the vector of the selected modes and SST
    DSM=PSM(:,1:N); %DSM is a binary upper triangular matrix
    MODES=PSM(:,end-1); %MODES is the proposed time duration between the 
                        %normal and the crash durations
    SST=PSM(:,end); %SST: Scheduled Start Time
    TD=zeros(N,1);
    CD=zeros(N,1);
    QD=zeros(N,1);
    RD=zeros(N,nR);
    for i=1:N
        if MODES(i)==0
            TD(i)=0;
            CD(i)=0;
            QD(i)=0;
            for j=1:nR
                RD(i,j)=0;
            end
        else
            TD(i)=T(i,MODES(i));
            CD(i)=C(i,MODES(i));
            QD(i)=q(i,MODES(i));
            for j=1:nR
                RD(i,j)=R(i,(j-1)*w+MODES(i));
            end
        end
    end
    tpts(I)=tptsst(DSM,TD,SST); %Calculate TPT for the member I of the 
                                %population
    tpcs(I)=diag(DSM)'*CD;      %Calculate TPC for the member I of the 
                                %population
    tpqs(I)=tpqfast(DSM,P,QD);  %Calculate TPQ for the member I of the 
                                %population
    c(I,1)=tpts(I)-Ct; %TPT should be lower than the time constraint
    c(I,2)=tpcs(I)-Cc; %TPC should be lower than the cost constraint
    c(I,3)=Cq-tpqs(I); %TPQ should be greater than the quality constraint
    c(I,4)=Cs-maxscore_PEM(DSM,P,Q); %TPS should be greater than Cs
    c(I,5:4+numel(CR))=maxresfun(SST,DSM,TD,RD)'-CR; %Resource constraint 
                                                     %should be satisfied
end
ceq=[];