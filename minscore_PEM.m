%Author: Zsolt T. Kosztyán Ph.D, University of Pannonia
%----------------
%Calculate minimal score value (PMAX) of possible project scenarios
%----------------
%Output:
%score: the minimal score value of the project scenario
%Inputs:
%PEM: N by N upper triangular adjacency matrix of logic network
%P: N by N score matrix of task/dependency inclusion
%Q: N by N score matrix of task/dependency exclusion
%----------------
%Usage:
%score=minscore_PEM(PEM,P,Q):
%----------------
%Example:
%
% %    |0.8, 0.4, 0.8|
% %PEM=|0.0, 0.7, 0.7| P=PEM Q=1-P
% %    |0.0, 0.0, 0.4|
%
%PEM=[[0.8,0.4,0.8];[0.0,0.7,0.7];[0.0,0.0,0.4]];P=PEM;Q=ones(3,3)-P;
%minscore_PEM(PEM,P,Q)

function score=minscore_PEM(PEM,P,Q)
score=1;
N=0;
p=diag(P);
q=diag(Q);
pem=diag(PEM);
N=numel(pem);   %HCs
pqmin=min([p,q],[],2);  %HCs
% for i=1:numel(pem)
%     if (p(i)>0)&&(p(i)<1)
%         N=N+1;
%         if pem(i)==1
%             score=score*p(i);
%         end
%         if pem(i)==0
%             score=score*q(i);
%         end
%         if (pem(i)<1)&&(pem(i)>0)
%             score=score*min(p(i),q(i));
%         end
%     end
% end
if N>0
%     score=power(score,1/N);
    score=power(prod([p(pem==1);q(pem==0);pqmin(pem>0&pem<1)]),1/N);    %HCs
end