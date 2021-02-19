%Author: Zsolt T. Kosztyán Ph.D, University of Pannonia
%----------------
%Calculate maximal score value (pMAX) for the set of  possible project 
%structures
%----------------
%Output:
%score: 1 by n vector of maximal score values of project structures
%Inputs:
%CELL_PEM: n element cell of N by N upper triangular adjacency matrices of 
%logic network
%P: N by N score matrix of task/dependency inclusion
%Q: N by N score matrix of task/dependency exclusion
%----------------
%Usage:
%score=core=maxscore_SNPM_cell(PEM,P,Q):
%----------------
%Example:
%
% %  |0.8, 0.4, 0.8|            |1 0 0|      |1 0 0|
% %P=|0.0, 0.7, 0.7| Q=1-P DSM1=|0 0 1| DSM2=|0 0 1|
% %  |0.0, 0.0, 0.4|            |0 0 0|      |0 0 1|
%
%P=[[1.0,0.4,0.8];[0.0,1.0,0.7];[0.0,0.0,1.0]];Q=ones(3,3)-P;
%DSM1=[[1,0,0];[0,0,1];[0,0,0]];DSM2=[[1,0,0];[0,0,1];[0,0,1]];
%maxscore_SNPM_cell({DSM1,DSM2},P,Q)

function score=maxscore_SNPM_cell(CELL_PEM,P,Q)
n=numel(CELL_PEM);
P=triu(P,1)+tril(P,-1); %Only out-diagonal values will be considered
Q=triu(Q,1)+tril(Q,-1); %Only out-diagonal values will be considered
score=zeros(1,n);
for i=1:n
    score(i)=maxscore_SNPM(CELL_PEM{i},P,Q);
end
