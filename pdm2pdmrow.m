%Author: Zsolt T. Kosztyán Ph.D, University of Pannonia
%----------------
%Convert an N by M PDM matrix to 1 by N*M row vector (for Excel users)
%----------------
%Output:
%pdmline: N*M row vector
%Input:
%PDM: N by M Project Domain Matrix 
%----------------
%Usage:
%pdmline=pdm2pdmrow(PDM)
%----------------
%Example:
%PDM=generatepdm(30,0.05,0,20,30,20,2,2,2);
%pdmline=pdm2pdmrow(PDM);

function pdmline=pdm2pdmrow(PDM)
pdmline=reshape(PDM,1,[]);