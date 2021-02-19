%Author: Zsolt T. Kosztyán Ph.D, University of Pannonia
%----------------
%Calculate maximal score values (PMAX) for the set of possible project 
%scenarios
%----------------
%Output:
%score: 1 by n vector of maximal score values of project scenarios
%Inputs:
%CELL_PEM: n element cell of N by N upper triangular adjacency matrices of 
%logic network
%P: N by N score matrix of task/dependency inclusion
%Q: N by N score matrix of task/dependency exclusion
%----------------
%Usage:
%score=maxscore_PEM_cell(CELL_PEM,P,Q)
%----------------
%Example:
%
% %  |0.8, 0.4, 0.8|            |1 0 0|      |1 0 0|
% %P=|0.0, 0.7, 0.7| Q=1-P DSM1=|0 0 1| DSM2=|0 0 1|
% %  |0.0, 0.0, 0.4|            |0 0 0|      |0 0 1|
%
%P=[[0.8,0.4,0.8];[0.0,0.7,0.7];[0.0,0.0,0.4]];Q=ones(3,3)-P;
%DSM1=[[1,0,0];[0,0,1];[0,0,0]];DSM2=[[1,0,0];[0,0,1];[0,0,1]];
%maxscore_PEM_cell({DSM1,DSM2},P,Q)

function score=maxscore_PEM_cell(CELL_PEM,P,Q)
n=numel(CELL_PEM);
score=zeros(1,n);
for i=1:n
    score(i)=maxscore_PEM(CELL_PEM{i},P,Q);
end