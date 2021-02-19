%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Calculate Cc=TPCmin+C%*(TPCmax-TPCmin)
%----------------
%Output:
%Cc=TPCmin+C%*(TPCmax-TPCmin)
%---------------- 
%Input:
%PDM=[PEM,TD,CD,{QD,RD}]: Project Domain Matrix
%where PEM is an N by N upper triangular matrix of logic domain, 
% TD is an N by w matrix of task durations
% CD is an N by w matrix of cost demands
% {QD is an N by w matrix of quality parameters} %optional
% {RD is an N by w*nR matrix of resource demands}%optional
%w=Number of completion modes
%ratio=C% cost ratio between interval [0,1]
%---------------- 
%Usage:
%Cc=percentc(PDM,w,ratio)
%---------------- 
%Example:
%PDM=generatepdm(30,0.05,0,20,30,20,2,2,2);
%w=2;
%ratio=0.5;
%Cc=percentc(PDM,w,ratio)

function Cc=percentc(PDM,w,ratio)
DSMdiag=ceil(diag(PDM(:,1:size(PDM,1)))); %All uncertain tasks/dependencies 
                                          %will be included
dsmdiag=floor(diag(PDM(:,1:size(PDM,1))));%All uncertain tasks/dependencies 
                                          %will be excluded
C=max(PDM(:,size(PDM,1)+w+1:size(PDM,1)+2*w),[],2);
c=min(PDM(:,size(PDM,1)+w+1:size(PDM,1)+2*w),[],2);
TPCmax=C'*DSMdiag;
TPCmin=c'*dsmdiag;
Cc=TPCmin+ratio*(TPCmax-TPCmin);