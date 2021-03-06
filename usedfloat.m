%Author: Zsolt T. Koszty�n Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Calculate total amount of used float
%----------------
%Output:
%UF is a scalar. The value is the total amount of used float
%---------------- 
%Input:
%SST: N by 1 vector of Scheduled Start Time
%---------------- 
%Usage:
%UF=usedfloat(SST)
%---------------- 
%Example:
%This is not executed as a stand-alone fuinction.
%---------------- 
%Prepositions and Requirements:
%1.) DSM matrix must be a binary upper triangular matrix
%2.) T should be a positive vector, R should be a positive matrix.

function UF=usedfloat(SST)
global DSM T
[~,EST]=tptfast(DSM,T);
SST=recalcsst(DSM,T,SST);
SST=reshape(SST,[],1);
UF=sum(SST-EST);