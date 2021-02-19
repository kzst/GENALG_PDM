function Multi_PDM=pdmtomultipdm(PDM,schedules)
N=numel(schedules);
n=size(PDM,1);
PEM=PDM(1:n,1:n);
m=size(PDM,2)-n;
Ds=PDM(:,n+1:n+m);
Multi_PDM=zeros((n+1)*N,(n+1)*N+m);
MN=(n+1)*N;
for i=1:N
    Multi_PDM((i-1)*(n+1)+1,(i-1)*(n+1)+1)=1; %New tasks
    Multi_PDM((i-1)*(n+1)+1,(i-1)*(n+1)+2)=1; %Connet PDM to the multiproject
    Multi_PDM(1,(i-1)*(n+1)+1)=1; %Connet PDM to the multiproject
   
    Multi_PDM((i-1)*(n+1)+2:i*(n+1),(i-1)*(n+1)+2:i*(n+1))=PEM;
    Multi_PDM((i-1)*(n+1)+2:i*(n+1),MN+1:end)=Ds;
    Multi_PDM((i-1)*(n+1)+1,MN+1:MN+2)=schedules(i);
end