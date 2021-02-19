%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Convert Project Schedule Matrix (PM) to Project Sctucture Matric (PSM) for
%RC-HDTQCTPs
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
%PM=[DSM,MODES,SST]: Project Schedule Matrix, where
% DSM is an N by N binary upper triangular matrix of logic domain
% MODES is an N by 1 vector selected completion modes
% SST is an N by 1 vector of scheduled starts
%TD is an N by w matrix of task durations
%CD is an N by w matrix of cost demands
%QD is an N by w matrix of quality paramters
%RD is an N by nR*w matrix of resource demands
%---------------- 
%Usage:
%PM=pmtopsmqr(PSM,TD,CD,QD,RD)
%---------------- 
%Example:
%This is not executed as a stand-alone fuinction. 

function PSM=pmtopsmqr(PM,TD,CD,QD,RD)
DSM=PM(:,1:size(PM,1));
w=size(TD,2);
nR=size(RD,2)/w;
MODES=PM(:,end-1);
N=numel(MODES);
t=zeros(N,1);
c=zeros(N,1);
q=zeros(N,1);
r=zeros(N,nR);
for i=1:N
    if MODES(i)==0
    else
        t(i)=TD(i,MODES(i));
        c(i)=CD(i,MODES(i));
        q(i)=QD(i,MODES(i));        
        for j=1:nR
            r(i,j)=RD(i,(j-1)*w+MODES(i));
        end
    end
end
PSM=[DSM,t,c,q,r,PM(:,end)];