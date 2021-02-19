%Author: Zsolt T. Koszty�n Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Generate a PDM matrix for TQCTP problems
%----------------
%Output:
%PDM=[PEM,TD,CD,QD,RD]: Project Domain Matrix
%where PEM is an N+nW by N+nW upper triangular matrix of logic domain, 
% TD is an N+nW by w matrix of task durations
% CD is an N+nW by w matrix of cost demands
% QD is an N+nW by w matrix of quality parameters
% RD is an N+nW by w*nR matrix of resource demands
%---------------- 
%Input:
%N=number of tasks
%ff=flexibility factor (between 0 to 1)
%cf=connectivity factor  (0,1,..)
%mT=max value of TD (must be positive)
%mC=max value of CD (must be positive)
%mR=max value of RD (must be positive)
%w=number of modes (1,2,..)
%nR=number of resources (1,2,..)
%nW=number of possible extra tasks (0,1,2,..)
%scale=Scale parameter
%---------------- 
%Usage:
%PDM=generatepdmq_scale(N,ff,cf,mTD,mCD,mRD,w,nR,nW)
%---------------- 
%Example:
%PDM=generatepdmq_scale(30,0.05,0,20,30,20,2,2,2,1.4)

function PDM=generatepdmq_scale(N,ff,cf,mTD,mCD,mRD,w,nR,nW,scale)
cf=cf+1;
PEM=phase3(triu(triu(min(ones(N)./max(repmat([1-cf:N-cf],N,1).^scale-...
    (repmat([0:N-1]',1,N).^scale),ones(N)),ones(N)),1)>rand(N),1)+...
    eye(N),ff,-.5); %Generate PEM matrix
nTD=w; %Width of TD = number of modes
nCD=w; %Width of CD = number of modes
nQD=w; %Width of CD = number of modes
nRD=w*nR; %Width of RD = number of modes x number of resources
TD=rand(N,nTD)*mTD; %Generate time domain
CD=rand(N,nCD)*mCD; %Generate cost domain
QD=rand(N,nQD);  %Generate quality domain
rD=rand(N,nRD)*mRD; %Generate resource domain
pem=zeros(N+nW);
pem(1:N,1:N)=PEM;
td=zeros(N+nW,nTD);
cd=zeros(N+nW,nCD);
qd=zeros(N+nW,nQD);
rd=zeros(N+nW,nRD);
if w==2 %In case of CTCTP the columns will be sorted
    TD=[max(TD,[],2)-max(TD,[],2).*rand(N,1)*.2,max(TD,[],2)];
    QD=[max(QD,[],2)-max(QD,[],2).*rand(N,1)*.2,max(QD,[],2)];    
    CD=[max(CD,[],2)-max(CD,[],2).*rand(N,1)*.2,max(CD,[],2)];
    RD=[];
    for i=1:2:nRD
        rmax=[max(rD(:,i:i+1,1)')'];
        rmin=rmax-[max(rD(:,i:i+1,1)')'].*rand(N,1)*.2;   
        RD=[RD,rmin,rmax];
    end
else
    RD=rD;
end
td(1:N,1:nTD)=TD;
cd(1:N,1:nCD)=CD;
qd(1:N,1:nQD)=QD;
rd(1:N,1:nRD)=RD;
PDM=[pem,td,cd,qd,rd];