%Author: Zsolt T. Koszty�n Ph.D, University of Pannonia
%----------------
%Calculate maximal score value (pMAX) of possible project structures
%----------------
%Output:
%score: the maximal score value of the project structure
%Inputs:
%PEM: N by N upper triangular adjacency matrix of logic network
%P: N by N score matrix of task/dependency inclusion
%Q: N by N score matrix of task/dependency exclusion
%----------------
%Usage:
%score=maxscore_SNPM(PEM,P,Q):
%----------------
%Example:
%
% %    |1.0, 0.4, 0.8|
% %PEM=|0.0, 1.0, 0.7| P=PEM Q=1-P
% %    |0.0, 0.0, 1.0|
%
%PEM=[[1.0,0.4,0.8];[0.0,1.0,0.7];[0.0,0.0,1.0]];P=PEM;Q=ones(3,3)-P;
%maxscore_SNPM(PEM,P,Q)

function score=maxscore_SNPM(PEM,P,Q)
score=1;
PEM=triu(PEM(:,1:size(PEM,1)));
N=0;
P=triu(P,1)+tril(P,-1); %Only out-diagonal values will be considered
Q=triu(Q,1)+tril(Q,-1); %Only out-diagonal values will be considered
SNPM=triu(PEM,1)+tril(PEM,-1);
for i=1:numel(SNPM(1,:))-1    
    for j=i:numel(SNPM(1,:))         
        if (P(i,j)>0)&&(P(i,j)<1)&&(i~=j)&&(PEM(i,i)>0)&&(PEM(j,j)>0)
            N=N+1; %Determine max score value
            if SNPM(i,j)==1 
                score=score*P(i,j);
            end
            if SNPM(i,j)==0
                score=score*Q(i,j);
            end
            if (SNPM(i,j)<1)&&(SNPM(i,j)>0)
                score=score*max(P(i,j),Q(i,j));
            end
        end
    end
end
if N>0 %The geometric mean of maximal out-diagonal score values.
    score=power(score,1/N);
end

