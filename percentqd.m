%Author: Zsolt T. Koszty�n Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Calculate Cq=TPQmin+Q%*(TPQmax-TPQmin)
%----------------
%Output:
%Cq=TPQmin+Q%*(TPQmax-TPQmin)
%---------------- 
%Input:
%PDM=[PEM,TD,CD,QD{,RD}]: Project Domain Matrix
%where PEM is an N by N upper triangular matrix of logic domain, 
% TD is an N by w matrix of task durations
% CD is an N by w matrix of cost demands
% QD is an N by w matrix of quality parameters %optional
% {RD is an N by w*nR matrix of resource demands}%optional
%w=Number of completion modes
%ratio=Q% quality ratio between interval [0,1]
%---------------- 
%Usage:
%Cq=percentq(PDM,w,ratio)
%---------------- 
%Example:
%PDM=generatepdmq(30,0.05,0,20,30,20,2,2,2);
%w=2;
%ratio=0.5;
%Cq=percentq(PDM,w,ratio)
%---------------- 
%Prepositions and Requirements:
%1.)Every score of inclusion should be greater than scores of exclusion

function Cq=percentqd(PDM,w,ratio)
N=size(PDM,1); %Number of tasks
PEM=PDM(:,1:N); %The original logic network
DSM=ceil(PEM); %If all uncertainties are realized
dsm=floor(PEM); %If all uncertainties are ignored
QD=PDM(:,N+2*w+1:N+3*w); % The quality domain
Q=max(QD,[],2); %The maximal quality level
q=min(QD,[],2); %The minimal quality level
TPQmax=tpqfastqd(DSM,PEM,QD,Q);
TPQmin=tpqfastqd(dsm,PEM,QD,q);
Cq=TPQmin+ratio*(TPQmax-TPQmin);