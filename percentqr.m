%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Calculate CR=TPRmin+R%*(TPRmax-TPRmin) when quality parameters are 
%considered
%----------------
%Output:
%CR=TPRmin+R%*(TPRmax-TPRmin)
%---------------- 
%Input:
%PDM=[PEM,TD,CD,QD,RD]: Project Domain Matrix
%where PEM is an N by N upper triangular matrix of logic domain, 
% TD is an N by w matrix of task durations
% CD is an N by w matrix of cost demands
% QD is an N by w matrix of quality parameters
% RD is an N by w*nR matrix of resource demands
%w=Number of completion modes
%ratio=R% ratio of resource demands between interval [0,1]
%---------------- 
%Usage:
%CR=percentqr(PDM,w,ratio)
%---------------- 
%Example:
%PDM=generatepdmq(30,0.05,0,20,30,20,2,2,2);
%w=2;
%ratio=0.5;
%CR=percentqr(PDM,w,ratio)'

function CR=percentqr(PDM,w,ratio)
DSM=floor(triu(PDM(:,1:size(PDM,1)),1))+diag(ceil(diag(PDM))); %If every 
   %tasks will be included, however, every dependencies will be excluded
dsm=ceil(triu(PDM(:,1:size(PDM,1)),1))+diag(floor(diag(PDM))); %If every 
   %tasks will be excluded, however, every dependencies will be included
rD=PDM(:,size(PDM,1)+3*w+1:end);
R=[]; %Maximal values of resource demands
r=[]; %Minimal values of resource demands
if w>1
    for i=1:w:size(rD,2)
        rmin=[min(rD(:,i:i+w-1,1),[],2)];
        rmax=[max(rD(:,i:i+w-1,1),[],2)];
        r=[r,rmin];
        R=[R,rmax];
    end
else
    R=rD;
    r=rD;
end
T=max(PDM(:,size(PDM,1)+1:size(PDM,1)+w),[],2); %min R when max T
t=min(PDM(:,size(PDM,1)+1:size(PDM,1)+w),[],2); %max R when min T
[~,EST,~,LST,~]=tptfast(DSM,t); %Optimization are within [EST,LST]
TPRmax=max(maxresfun(EST,DSM,t,R),maxresfun(LST,DSM,t,R));
if ratio==1.0
    CR=TPRmax;
else
    SST=paretores(dsm,T,r); %calculationof TPRmin
    TPRmin=maxresfun(SST,dsm,T,r);
    CR=TPRmin+ratio*(TPRmax-TPRmin);
end