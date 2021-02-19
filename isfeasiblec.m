%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Feasibility check for crossover function of GA approach for HCTCTPs
%----------------
%Output:
%loout=logical value of the feasibility
%---------------- 
%Inputs:
%DSM=N by N upper trisngle binary matrix of logic domain
%MODES=N by 1 vector of completion modes (i.e. duration time)
%---------------- 
%Usage:
%loout=isfeasiblec(DSM,MODES)
%---------------- 
%Example:
%This is not executed as a stand-alone fuinction. Instead of specifying
%global values this function is cannot be run.
%---------------- 
%Prepositions and Requirements:
%1.)Global values must be specified.

function loout=isfeasiblec(DSM,MODES)
global P T C Ct Cc Cs
%global values: 
% P is the PopSize by PopSize scores matrix of task/dependency inclusion
% T is the column vector of time demands (=the time domain)
% C is the column vector of cost demands (=the cost domain)
% Ct is the time constraint
% Cc is the cost constraint
% Cs is the score constraint

loout=false;
TD=zeros(numel(MODES),1); %Calculation of realized time demands
CD=zeros(numel(MODES),1); %Calculation of realized cost demands
for i=1:numel(MODES)
    if DSM(i,i)==0
    else
        TD(i)=MODES(i);
        CD(i)=timetocost(TD(i),T(i,:),C(i,:));
    end
end
TPS=maxscore_PEM(DSM,P,ones(size(P,1))-P);
TPT=tptfast(DSM,TD);
TPC=tpcfast(DSM,CD);
if max([TPT-Ct,TPC-Cc,Cs-TPS])<=0 %The TPT should be lower than Ct, TPC 
    loout=true; %should be lower than Cc and TPS should be greater than Cs
end