%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Calculate the set of total amount of used floats (GPU READY VERSION)
%----------------
%Output:
%UF is an 1 by n vector of total amounts of used floats
%---------------- 
%Input:
%CELL_PM: n element of project schedule matrix (PM=[DSM,MODES,SST]), where
%DSM: N by N binary upper triangular matrix of the logic domain
%MODES is an N by 1 vector of completion modes
%SST: N by 1 vector of Scheduled Start Time
%---------------- 
%Usage:
%UF=usedfloatfast_cell(CELL_PM)
%---------------- 
%Example:
%This is not executed as a stand-alone fuinction.
%---------------- 
%Prepositions and Requirements:
%1.) T should be a positive vector, R should be a positive matrix.

function UF=usedfloatfast_cell(CELL_PM)
global T
n=numel(CELL_PM);
w=size(T,2);
UF=zeros(n,1);
for i=1:n
    dsm=CELL_PM{i};
    N=size(dsm,1);
    MODES=dsm(:,end-1);
    SST=dsm(:,end);
    t=zeros(N,1);
    for I=1:w
        if MODES(I)==0
            t(I)=0;            
        else
            t(I)=T(I,MODES(I));           
        end
    end
    dsm=dsm(:,1:N);
    [~,EST]=tptfast(dsm,t);
    SST=recalcsstfast(dsm,t,SST);
    SST=reshape(SST,[],1);
    UF(i)=sum(SST-EST);
end