%Author: Zsolt T. Kosztyan Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Calculate maximum resource demands
%----------------
%Output:
%rMAX: is an nR by 1 vector of maximum resource demands
%---------------- 
%Inputs:
%SST: N by 1 vector of Scheduled Start Time
%---------------- 
%Usage:
%rMAX=maxres(SST)
%---------------- 
%Example:
%This is not executed as a stand-alone fuinction. Instead of specifying
%global values this function is cannot be run.
%---------------- 
%Prepositions and Requirements:
%1.) DSM matrix must be a binary upper triangular matrix
%2.) SST,T should be a positive vector, R should be a positive matrix.

function rMAX=maxres(SST)
global DSM T R
rMAX=maxresfun(SST,DSM,T,R);