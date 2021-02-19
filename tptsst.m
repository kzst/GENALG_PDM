%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Evaluate the project duration, when the Scheduled Start Time (SST) is 
%specified
%----------------
%Output:
%TPT: Total Project Time (scalar)
%SST: Scheduled Start Time (N by 1 vector)
%SFT: Scheduled Finish Time (N by 1 vector)
%----------------
%Inputs:
%DSM: Dependeny Structure Matrix (N by N matrix (logic plan))
%T: Duration time (N by 1 vector)
%SST: Scheduled Start Time (N by 1 vector)
%----------------
%Usage: 
%[TPT,SST,SFT]=tptsst(DSM,T,SST)
%----------------
%Example: - Comparison of TPT when task scheduled at as early / as late as 
%possible
% %    |1,1,0|    |3|
% %DSM=|0,1,0|, T=|4|
% %    |0,0,1|    |5|
%
%DSM=[[1,1,0];[0,1,0];[0,0,1]]; %Specification of DSM
%T=[3,4,5]';
%[TPT,SST,SFT]=tptsst(DSM,T,[1,1,1]')
function [TPT,SST,SFT]=tptsst(DSM,T,SST)
SST=reshape(SST,[],1);
N=numel(T); %Number of elements in vector T
TPT=0;
if numel(find(diag(DSM)>0))
    %Initialisation
    SFT=SST+T; %EFT=EST+T
    %Forward pass
    for i=1:(N-1)
        for j=(i+1):N
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
    TPT=max(SFT);
end