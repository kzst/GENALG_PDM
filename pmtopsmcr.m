%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Convert Project Schedule Matrix (PM) to Project Sctucture Matric (PSM) for
%RC-HCTCTPs
%----------------
%Output:
%PSM=[DSM,t,c,q,r,SST], where
% DSM is an N by N binary upper triangular matrix of logic domain
% t is an N by 1 vector of time demands
% c is an N by 1 vector of cost demands
% r is an N by nR matrix of resource parameters
% SST is an N by 1 vector of scheduled starts
%---------------- 
%Inputs:
%PM=[DSM,T,SST]: Project Schedule Matrix, where
% DSM is an N by N binary upper triangular matrix of logic domain
% T is an N by 1 vector of realized task durations
% SST is an N by 1 vector of scheduled starts
%TD is an N by w matrix of task durations
%CD is an N by w matrix of cost demands
%RD is an N by nR*w matrix of resource demands
%---------------- 
%Usage:
%PM=pmtopsmcr(PSM,TD,CD,RD)
%---------------- 
%Example:
%This is not executed as a stand-alone fuinction. 

function PSM=pmtopsmcr(PM,TD,CD,RD)
DSM=PM(:,1:size(PM,1));
w=size(TD,2);
nR=size(RD,2)/w;
MODES=PM(:,end-1);
N=numel(MODES);
t=zeros(N,1);
c=zeros(N,1);
r=zeros(N,nR);
for i=1:N
    if MODES(i)==0
    else
        t(i)=MODES(i);
        c(i)=timetocost(t(i),TD(i,:),CD(i,:));
        r(i,:)=timetor(t(i),TD(i,:),RD(i,:));
    end
end
PSM=[DSM,t,c,r,PM(:,end)];