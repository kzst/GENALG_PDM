%Author: Zsolt T. Kosztyán Ph.D, University of Pannonia
%----------------
%Convert an N*M row matrix to a PDM matrix (for Excel users)
%----------------
%Output:
%PDM: N by M Project Domain Matrix 
%Input:
%pdmline: N*M row vector
%----------------
%Usage:
%PDM=pdmrow2pdm(pdmline,n,m)
%----------------
%Example:
%PDM=generatepdm(30,0.05,0,20,30,20,2,2,2);
%pdmline=pdm2pdmrow(PDM); %PDM=>pdmline
%PDM=pdmrow2pdm(pdmline,size(PDM,1),size(PDM,2)); %pdmline=>PDM

function PDM=pdmrow2pdm(pdmline,n,m)
PDM=reshape(pdmline,n,m);