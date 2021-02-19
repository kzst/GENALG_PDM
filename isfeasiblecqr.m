%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Feasibility check for crossover function of GA approach for 
%resource-constrained HCTQCTPs
%----------------
%Output:
%loout=logical value of the feasibility
%---------------- 
%Inputs:
%DSM=N by N upper trisngle binary matrix of logic domain
%MODES=N by 1 vector of completion modes
%SST=N by 1 vector of scheduled strart time
%---------------- 
%Usage:
%loout=isfeasiblecqr(DSM,MODES,SST)
%---------------- 
%Example:
%This is not executed as a stand-alone fuinction. Instead of specifying
%global values this function is cannot be run.
%---------------- 
%Prepositions and Requirements:
%1.)Global values must be specified.

function loout=isfeasiblecqr(DSM,MODES,SST)
global P T C q R Ct Cc Cq CR Cs
%global values: 
% P is the PopSize by PopSize scores matrix of task/dependency inclusion
% T is the column vector of time demands (=the time domain)
% C is the column vector of cost demands (=the cost domain)
% R is the N by nR matrix of resource constraint (=the resource domain)
% q is the column vector of quality parameters (=the quality domain)
% Ct is the time constraint
% Cc is the cost constraint
% CR is the 1 by nR row vector of resource constraint
% Cq is the quality constraint

loout=false;
TD=zeros(numel(MODES),1); %Calculation of realized time demands
CD=zeros(numel(MODES),1); %Calculation of realized cost demands
QD=zeros(numel(MODES),1); %Calculation of realized quality demands
RD=zeros(numel(MODES),numel(CR)); %Calculation of realized resource demands
for i=1:numel(MODES)
    if DSM(i,i)==0
    else
        TD(i)=MODES(i);
        CD(i)=timetocost(TD(i),T(i,:),C(i,:));
        QD(i)=timetoquality(TD(i),T(i,:),q(i,:));        
        RD(i,:)=timetor(TD(i),T(i,:),R(i,:));
    end
end
TPS=maxscore_PEM(DSM,P,ones(size(P,1))-P);
TPT=tptsst(DSM,TD,SST);
TPC=tpcfast(DSM,CD);
TPQ=tpqfast(DSM,P,QD);
TPR=maxresfun(SST,DSM,TD,RD)';
if max([TPT-Ct,TPC-Cc,Cq-TPQ,TPR-CR,Cs-TPS])<=0 %The TPT should be lower   
    loout=true; %than Ct, TPC should be lower than Cc and TPS should be  
end %greater than Cs and TPQ should be greater than Cq and TPR for every 
%resources should be lower than CR 