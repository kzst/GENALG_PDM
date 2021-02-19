%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Convert Project Schedule Matrix (PM) to Project Sctucture Matric (PSM) for
%HCTCTPs
%----------------
%Output:
%PSM=[DSM,t,c,EST], where
% DSM is an N by N binary upper triangular matrix of logic domain
% t is an N by 1 vector of time demands
% c is an N by 1 vector of cost demands
% EST is an N by 1 vector of earliest schedule
%---------------- 
%Inputs:
%PM=[DSM,T]: Project Schedule Matrix, where
% DSM is an N by N binary upper triangular matrix of logic domain
% T is an N by 1 vector of realized task durations
%TD is an N by w matrix of task durations
%CD is an N by w matrix of cost demands
%---------------- 
%Usage:
%PSM=pmtopsmc(PM,TD,CD)
%---------------- 
%Example:
%This is not executed as a stand-alone fuinction

function PSM=pmtopsmc(PM,TD,CD)
DSM=PM(:,1:size(PM,1));
MODES=PM(:,end);
N=numel(MODES);
t=zeros(N,1);
c=zeros(N,1);
for i=1:N
    if MODES(i)==0
    else
        t(i)=MODES(i);
        c(i)=timetocost(t(i),TD(i,:),CD(i,:));
    end
end
[~,EST]=tptfast(DSM,t);
PSM=[DSM,t,c,EST];