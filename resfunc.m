%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Calculate breakpoints and changes in resource demands
%----------------
%Outputs:
%BP is a row vector of break points
%RESFUNC is a row vector of resource values 
%---------------- 
%Inputs:
%DSM: N by N binary upper triangular matrix of the logic domain
%SST: N by 1 vector of Scheduled Start Time
%T is an N by 1 vector of task durations
%R is an N by nR matrix of resource demands
%---------------- 
%Usage:
%[BP,RESFUNC]=resfunc(DSM,SST,T,R)
%---------------- 
%Example: draw resource demands for tasks, which scheduled as early as
%possible
%DSM=triu(round(rand(10))); T=rand(10,1)*20; R=rand(10,3)*5; %generate
%initial domains
%[~,EST,~,~,~]=tptfast(DSM,T); %calculate EST
%[BP,RESFUNC]=resfunc(DSM,EST,T,R); %calculate BP and RESFUNC
%figure('Name','EST');
%drawres(BP,RESFUNC);
%---------------- 
%Prepositions and Requirements:
%1.) DSM matrix must be a binary upper triangular matrix
%2.) T should be a positive vector, R should be a positive matrix.

function [BP,RESFUNC]=resfunc(DSM,SST,T,R)
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
BP=sort(union(SST,SFT));
b=numel(BP); %number of breakpoints
r=numel(R(1,:)); %number of resources
RESFUNC=zeros(b,r);
for i=1:b
    RESFUNC(i,:)=sum(R(find((SST<=BP(i))&(SFT>BP(i))),:),1); %calculate resource function
end
