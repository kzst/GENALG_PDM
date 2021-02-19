function pdmtopdmq(INFILE,OUTFILE)
load(INFILE);
N=size(PDMs,1);
M=size(PDMs,2);
K=size(PDMs,3);
C=size(sCONSs,2);
PDMsq=zeros(N,M+w,K);
sSCONSsq=zeros(K,C+1);
cSCONSsq=zeros(K,C+1);
for i=1:K
    i
    Q=rand(N,w);
    if w==2
        Q=[min(Q,[],2),max(Q,[],2)];
    end
    PDM=[PDMs(:,1:N+2*w,i),Q,PDMs(:,N+2*w+1:end,i)];
    PDMsq(:,:,i)=PDM;
    sCpq=PROPs(i,7);
    cCpq=PROPs(i,11);
    sSCONSq=[sCONSs(i,1:2),percentq(PDM,w,sCpq),sCONSs(i,3:end)];
    cSCONSq=[cCONSs(i,1:2),percentq(PDM,w,cCpq),cCONSs(i,3:end)];
    sSCONSsq(i,:)=sSCONSq;
    cSCONSsq(i,:)=cSCONSq;
end
PDMs=PDMsq;
sCONSs=sSCONSsq;
cCONSs=cSCONSsq;
save(OUTFILE,'cCONSs','cf','mCD','mRD','mTD','nR','PDMs','PROPs',...
    'sCONSs','w')