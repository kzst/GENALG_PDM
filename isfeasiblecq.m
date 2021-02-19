%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Feasibility check for crossover function of GA approach for HCTQCTPs
%----------------
%Output:
%loout=logical value of the feasibility
%---------------- 
%Inputs:
%DSM=N by N upper trisngle binary matrix of logic domain
%MODES=N by 1 vector of completion modes
%---------------- 
%Usage:
%loout=isfeasiblecq(DSM,MODES)
%---------------- 
%Example:
%This is not executed as a stand-alone fuinction. Instead of specifying
%global values this function is cannot be run.
%---------------- 
%Prepositions and Requirements:
%1.)Global values must be specified.

function loout=isfeasiblecq(DSM,MODES)
global P T C q Ct Cc Cq Cs
%global values: 
% P is the PopSize by PopSize scores matrix of task/dependency inclusion
% T is the column vector of time demands (=the time domain)
% C is the column vector of cost demands (=the cost domain)
% q is the column vector of quality parameters (=the quality domain)
% Ct is the time constraint
% Cc is the cost constraint
% Cs is the score constraint
% Cq is the quality constraint

loout=false;
TD=zeros(numel(MODES),1); %Calculation of realized time demands
CD=zeros(numel(MODES),1); %Calculation of realized cost demands
QD=zeros(numel(MODES),1); %Calculation of realized quality demands
for i=1:numel(MODES)
    if DSM(i,i)==0
    else
        TD(i)=MODES(i);
        CD(i)=timetocost(TD(i),T(i,:),C(i,:));
        QD(i)=timetoquality(TD(i),T(i,:),q(i,:));        
    end
end
TPS=maxscore_PEM(DSM,P,ones(size(P,1))-P);
TPT=tptfast(DSM,TD);
TPC=tpcfast(DSM,CD);
TPQ=tpqfast(DSM,P,QD);
if max([TPT-Ct,TPC-Cc,Cq-TPQ,Cs-TPS])<=0%The TPT should be lower than Ct,  
    loout=true; %TPC should be lower than Cc and TPS should be greater 
end %than Cs and TPQ should be greater than Cq