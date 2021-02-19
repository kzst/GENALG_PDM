%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%GPU ready version of recalculation SST
%----------------
%Output:
%SST N by 1 vector of scheduled start time
%---------------- 
%Inputs:
%DSM N by N uppert triangular binary matrix of logic domai
%T N by 1 vector of task durations
%SST N by 1 vector of initial scheduled start time
%---------------- 
%Usage:
%SST=recalcsstfast(DSM,T,SST))
%---------------- 
%Example:
%This is not executed as a stand-alone fuinction.

function SST=recalcsstfast(DSM,T,SST)
n=numel(T); %Number of elements in vector T
SST=reshape(SST,[],1);
n=numel(T); %Number of elements in vector T
T(diag(DSM)==0)=0;
SFT=SST+T;
dsm=triu(DSM,1);
dsm(diag(DSM)==0,:)=0; 
dsm(:,diag(DSM)==0)=0;
for i=1:n
    SST=max(bsxfun(@times,dsm,SFT))';
    SFT=SST+T;
end