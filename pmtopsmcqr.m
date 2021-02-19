%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Convert Project Schedule Matrix (PM) to Project Sctucture Matric (PSM) for
%RC-HCTQCTPs
%----------------
%Output:
%PSM=[DSM,t,c,q,r,SST], where
% DSM is an N by N binary upper triangular matrix of logic domain
% t is an N by 1 vector of time demands
% c is an N by 1 vector of cost demands
% q is an N by 1 vector of quality paramters
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
%QD is an N by w matrix of quality paramters
%RD is an N by nR*w matrix of resource demands
%---------------- 
%Usage:
%PM=pmtopsmcqr(PSM,TD,CD,QD,RD)
%---------------- 
%Example:
%This is not executed as a stand-alone fuinction. 

function PSM=pmtopsmcqr(PM,TD,CD,QD,RD)
DSM=PM(:,1:size(PM,1));
M=size(TD,2);
nR=size(RD,2)/M;
MODES=PM(:,end-1);
n=numel(MODES);
t=zeros(n,1);
c=zeros(n,1);
q=zeros(n,1);
r=zeros(n,nR);
for i=1:n
    if MODES(i)==0
    else
        t(i)=MODES(i);
        c(i)=timetocost(t(i),TD(i,:),CD(i,:));
        q(i)=timetoquality(t(i),TD(i,:),QD(i,:));
        r(i,:)=timetor(t(i),TD(i,:),RD(i,:));
    end
end
PSM=[DSM,t,c,q,r,PM(:,end)];