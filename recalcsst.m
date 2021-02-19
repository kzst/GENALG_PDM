%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Recalculate SST based on a given DSM
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
%SST=recalcsst(DSM,T,SST))
%---------------- 
%Example:
%This is not executed as a stand-alone fuinction.

function SST=recalcsst(DSM,T,SST)
n=numel(T); %Number of elements in vector T
SST=reshape(SST,[],1);
if numel(find(diag(DSM)>0))
    %Initialisation
    SFT=SST+T; %EFT=EST+T
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
end