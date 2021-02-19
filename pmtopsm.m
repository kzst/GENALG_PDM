%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Convert Project Schedule Matrix (PM) to Project Sctucture Matric (PSM) for
%HDTCTPs
%----------------
%Output:
%PSM=[DSM,t,c,EST], where
% DSM is an N by N binary upper triangular matrix of logic domain
% t is an N by 1 vector of time demands
% c is an N by 1 vector of cost demands
% EST is an N by 1 vector of earliest schedule
%---------------- 
%Inputs:
%PM=[DSM,MODES]: Project Schedule Matrix, where
% DSM is an N by N binary upper triangular matrix of logic domain
% MODES is an N by 1 vector of selected completion modes
%TD is an N by w matrix of task durations
%CD is an N by w matrix of cost demands
%---------------- 
%Usage:
%PM=pmtopsm(PSM,TD,CD)
%---------------- 
%Example:
%This is not executed as a stand-alone fuinction.

function PSM=pmtopsm(PM,TD,CD)
DSM=PM(:,1:size(PM,1));
MODES=PM(:,end);
N=numel(MODES);
t=zeros(N,1);
c=zeros(N,1);
for i=1:N
    if MODES(i)==0
    else
        t(i)=TD(i,MODES(i));
        c(i)=CD(i,MODES(i));        
    end
end
[~,EST]=tptfast(DSM,t);
PSM=[DSM,t,c,EST];