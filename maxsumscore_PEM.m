%Author: Zsolt T. Kosztyán Ph.D, University of Pannonia
%----------------
%Calculate maximal sum of score values of possible project scenarios
%----------------
%Output:
%score: the maximal score value of the project scenario
%Inputs:
%PEM: N by N upper triangular adjacency matrix of logic network
%P: N by N score matrix of task/dependency inclusion
%Q: N by N score matrix of task/dependency exclusion
%----------------
%Usage:
%maxsumscore_PEM(PEM,P,Q)
%----------------
%Example:
%
% %    |0.8, 0.4, 0.8|
% %PEM=|0.0, 0.7, 0.7| P=PEM Q=1-P
% %    |0.0, 0.0, 0.4|
%
%PEM=[[0.8,0.4,0.8];[0.0,0.7,0.7];[0.0,0.0,0.4]];P=PEM;Q=ones(3,3)-P;
%maxsumscore_PEM(PEM,P,Q)

function score=maxsumscore_PEM(PEM,P,Q)
score=1;
p=diag(P);
q=diag(Q);
pem=diag(PEM);
N=numel(pem);   
pqmax=max([p,q],[],2);  
if N>0 %The score of the project scenario is the sum of maximum task scores
    score=sum([p(pem==1);q(pem==0);pqmax(pem>0&pem<1)]); 
end                                                               