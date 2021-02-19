%Author: Zsolt T. Kosztyan Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Feasibility check of PSMs
%----------------
%Output:
%L   : Logical value (1 if feasible, 0 if infeasible)
%TPT : Total Project Time (scalar)
%TPC : Total Project Cost (scalar)
%TPR : Total Project Resources (vector)
%TPS : Total Project Score/Scope (scalar)
%PSM: N by M+1 matrix of the calculated project plan. PSM contains the
%----------------
%Inputs:
%PDM: N by N+(M-N)*W matrix of the stochastic project plan. 
%PSM=[DSM,TD,CD{,RD},EST|SST] matrix of output structure (Project Structure
%Matrix)


function [L,TPT,TPC,TPR,TPS]=feasibilitycheck(PSM,PDM,c)
L=true;
DSM=PSM(:,1:size(PSM));
PEM=PDM(:,1:size(PSM));
p=PEM;
q=ones(size(PEM))-p;
T=PSM(:,size(PSM,1)+1);
C=PSM(:,size(PSM,1)+2);
SST=PSM(:,end);
TPT=tptfast(DSM,T);
TPC=tpcfast(DSM,C);
TPS=maxscore_PEM(DSM,p,q);
if size(PSM,1)+3>size(PSM,2) %If it has no resources
    R=PSM(:,size(PSM,1)+3:end-1);
    TPR=maxresfun(SST,DSM,T,R)';
    L=min([TPT,TPC,TPR,-TPS]<=[c(1),c(2),c(3:end-1),c(end)]);
else
    R=[];
    TPR=[];
    L=min([TPT,TPC,-TPS]<=[c(1),c(2),c(end)]);
end

