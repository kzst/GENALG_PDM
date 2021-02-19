%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Calculate Cs=TPSmin+S%*(TPSmax-TPSmin)
%----------------
%Output:
%Cs=TPSmin+S%*(TPSmax-TPSmin)
%---------------- 
%Input:
%PDM=[PEM,TD,CD,{QD,RD}]: Project Domain Matrix
%where PEM is an N by N upper triangular matrix of logic domain, 
% TD is an N by w matrix of task durations
% CD is an N by w matrix of cost demands
% {QD is an N by w matrix of quality parameters} %optional
% {RD is an N by w*nR matrix of resource demands}%optional
%ratio=S% score ratio of inclusion, this value should be between interval 
%[0,1]
%---------------- 
%Usage:
%Cs=percents(PDM,ratio)
%---------------- 
%Example:
%PDM=generatepdm(30,0.05,0,20,30,20,2,2,2);
%ratio=0.5;
%Cs=percents(PDM,ratio)

function Cs=percents(PDM,ratio)
PEM=PDM(:,1:size(PDM,1)); %N by N matrix of the logic domain
TPSmax=maxscore_PEM(PEM,PEM,ones(size(PEM,1))-PEM);
TPSmin=minscore_PEM(PEM,PEM,ones(size(PEM,1))-PEM);
Cs=TPSmin+ratio*(TPSmax-TPSmin);