%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Calculate Ct=TPTmin+T%*(TPTmax-TPTmin)
%----------------
%Output:
%Ct=TPTmin+T%*(TPTmax-TPTmin)
%---------------- 
%Input:
%PDM=[PEM,TD,CD,{QD,RD}]: Project Domain Matrix
%where PEM is an N by N upper triangular matrix of logic domain, 
% TD is an N by w matrix of task durations
% CD is an N by w matrix of cost demands
% {QD is an N by w matrix of quality parameters} %optional
% {RD is an N by w*nR matrix of resource demands}%optional
%w=Number of completion modes
%ratio=T% task duration ratio between interval [0,1]
%---------------- 
%Usage:
%Ct=percentt(PDM,w,ratio)
%---------------- 
%Example:
%PDM=generatepdm(30,0.05,0,20,30,20,2,2,2);
%w=2;
%ratio=0.5;
%Ct=percentt(PDM,w,ratio)

function Ct=percentt(PDM,w,ratio)
DSM=ceil(PDM(:,1:size(PDM,1))); %All uncertain tasks/dependencies 
                                      %will be included
dsm=floor(PDM(:,1:size(PDM,1)));%All uncertain tasks/dependencies 
                                      %will be excluded
T=max(PDM(:,size(PDM,1)+1:size(PDM,1)+w),[],2);
t=min(PDM(:,size(PDM,1)+1:size(PDM,1)+w),[],2);
TPTmax=tptfast(DSM,T);
TPTmin=tptfast(dsm,t);
Ct=TPTmin+ratio*(TPTmax-TPTmin);