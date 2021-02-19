%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Calculate average resource demands
%----------------
%Output:
%rAVG: is an nR by 1 vector of average resource demands
%---------------- 
%Inputs:
%SST: N by 1 vector of Scheduled Start Time
%DSM: N by N binary upper triangular matrix of the logic domain
%T is an N by 1 vector of task durations
%R is an N by nR matrix of resource demands
%---------------- 
%Usage:
%rAVG=avgresfun(SST,DSM,T,R)
%---------------- 
%Example:
%DSM=round(triu(rand(10))); %Generate logic domain
%T=rand(10,1)*20; %Generate logic domain
%[~,EST]=tptfast(DSM,T); %Calculate Early Start Times
%R=rand(10,3)*5; %Generate resource demands (nR=3)
%rAVG=avgresfun(EST,DSM,T,R) %Clculate average resource demands for the
%       %three resources
%---------------- 
%Prepositions and Requirements:
%1.) DSM matrix must be a binary upper triangular matrix
%2.) SST,T should be a positive vector, R should be a positive matrix.
function rAVG=avgresfun(SST,DSM,T,R)
N=numel(SST); %Number of tasks
DSM=round(triu(DSM)); %DSM must be an upper triangular binary matrix.
DSM(diag(DSM)==0,:)=0; %Excluded task has no dependency
DSM(:,diag(DSM)==0)=0; %Excluded task has no dependency
T(diag(DSM)==0)=0;  %Excluded task has no time demands
R(diag(DSM)==0,:)=0; %Excluded task has no resources
R(isnan(R))=0;
SST(diag(DSM)==0)=0; %Excluded task's start time is zero
SST=reshape(SST,N,1); %%SST should be a column vector
n=numel(T); %Number of elements in vector T
%Initialisation
SFT=SST+T; 
%Forward pass
for i=1:(n-1)
    for j=(i+1):n
        if DSM(i,i)>0
            if DSM(j,j)>0
                if DSM(i,j)>0 %If there is a dependency between task i and task j...
                    if SST(j)<SFT(i)
                        SST(j)=SFT(i);
                        SFT(j)=SST(j)+T(j);
                    end
                end
            end
        end
    end
end
BP=sort(union(SST,SFT)); %Breakpoints, where the resource demands should 
%be recalculated
b=numel(BP); %number of breakpoints
RESFUNC=mtimes(((repmat(SST,1,b)<=repmat(BP',N,1))&(repmat(SFT,...
    1,b)>repmat(BP',N,1)))',R); %Calculate resource function
rAVG=mean(RESFUNC)'; %Calculkate average of resource demands