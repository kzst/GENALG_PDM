%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Evaluate ratio of flexible dependency
%----------------
%Output:
%Fr = Ratio of flexible dependency
%----------------
%Input:
%PDM = Project Domain Matrix
%----------------
%Usage: 
%Fr=fratio(PDM)
%----------------
%Example: - Calculate flexible dependencies
%
%PDM=generatepdm(30,0.05,0,20,30,20,2,2,2)
%Fr=fratio(PDM)

function Fr=fratio(PDM)
pem=triu(PDM(:,1:size(PDM,1)),1);
if numel(pem(pem>0))>0
    Fr=numel(pem(pem>0&pem<1))/numel(pem(pem>0));
else
    Fr=-1;
end