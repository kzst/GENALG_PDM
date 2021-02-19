%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Calculate cost demands of a project scenario
%----------------
%Output:
%TPC: Total Project Cost (scalar)
%---------------- 
%Inputs:
%DSM: Upper triangular binary matrix of logic domain
%C:   N by 1 vector of cost demands
%---------------- 
%Usage:
%TPC=tpcfast(DSM,C)
%---------------- 
%Example: Calculate TPC for a generated project scenario
%DSM=triu(round(rand(10)*.5+.5)); %Genarate DSM
%C=rand(10,1)*30; %Generate C vector (cost demands)
%TPC=tpcfast(DSM,C) %Calculate TPC
%---------------- 
%Prepositions and Requirements:
%1.)The diagonal of DSM matrix must be binary. (=DSM is a project scenario)
%2.)C contains real numbers

function TPC=tpcfast(DSM,C)
TPC=diag(ceil(DSM))'*C; %The diagonal of DSM must be a binary matrix, if 
%not, we consider the round towards plus infinity of the DSM's diagonal