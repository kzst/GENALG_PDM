%Author: Zsolt T. Koszty�n Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Evaluate EST, EFT, LST and LFT times of activity
%Evaluate the project duration. This code contains only 1-depth for cycles
%in order to help the parallelization
%----------------
%Outputs:
%TPT: Total Project Time (scalar)
%EST: Early Start Time (N by 1 vector)
%EFT: Early Finish Time (N by 1 vector)
%LST: Latest Start Time (N by 1 vector)
%LFT: Latest Finish Time (N by 1 vector)
%----------------
%Inputs:
%DSM: Dependeny Structure Matrix (N by N matrix (logic plan))
%T: Duration time (N by 1 vector)
%----------------
%Usage: 
%[TPT,EST,EFT,LST,LFT]=tptfast(DSM,T)
%----------------
%Example: - Calculate TPT and EST,LST,EFT,LFTs
% %    |1,1,0|    |3|
% %DSM=|0,1,0|, T=|4|
% %    |0,0,1|    |5|
%
%[TPT,EST,EFT,LST,LFT]=tptfast([[1,1,0];[0,1,0];[0,0,1]],[3,4,5]')

function [TPT,EST,EFT,LST,LFT]=tptfast(DSM,T)
DSM=round(DSM); %DSM must be binary matrix
T=real(T); %T must be real
N=numel(T); %Number of elements in vector T
EST=zeros(N,1); %EST will be an N by 1 null vector at the first step
T(diag(DSM)==0)=0; %The escluded task duration is irrelevant=>set to be 0
T=reshape(T,[],1); %T must be column vector
EFT=EST+T; %EFTi=ESTi+Ti (i=1..N)
dsm=triu(DSM,1);
dsm(diag(DSM)==0,:)=0; 
dsm(:,diag(DSM)==0)=0;
for i=1:N %Forward pass
    EST=max(bsxfun(@times,dsm,EFT))'; %EST calculated step by step
    EFT=EST+T; %EFTi=ESTi+Ti (i=1..N)
end
TPT=max(EFT); %TPT is the makespan of the longest path
LFT=repmat(TPT,N,1); %LFTi=TPT (i=1..N)
LST=LFT-T; %LSTi=LFTi-Ti (i=1..N)
for i=1:N %Backward pass
    TPTDSM=bsxfun(@times,dsm',LST)'; %Calculate LST step by step
    TPTDSM(isnan(TPTDSM))=TPT; %Independent tasks' LFT = TPT
    TPTDSM(~TPTDSM)=TPT; %Independent tasks' LFT = TPT
    LFT=min(TPTDSM,[],2); %Calculate LFT
    LST=LFT-T; %LSTi=LFTi-Ti (i:=1..N)
end