%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Evaluate ratio of supplementary tasks
%----------------
%Output:
%Sr = Ratio of supplementary tasks
%----------------
%Input:
%PDM = Project Domain Matrix
%----------------
%Usage: 
%Sr=sratio(PDM)
%----------------
%Example: - Calculate flexible dependencies
%
%PDM=generatepdm(30,0.05,0,20,30,20,2,2,2)
%Sr=sratio(PDM)

function Sr=sratio(PDM)
pem=diag(PDM(:,1:size(PDM,1)));
if numel(pem(pem>0))>0
    Sr=numel(pem(pem>0&pem<1))/numel(pem(pem>0));
else
    Sr=-1;
end