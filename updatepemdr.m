%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Subfunction for conversion between chromosome and PSM
%----------------
%Output:
%PSM: An N by N+2 PSM=[DSM,MODES,SST] matrix, where
% DSM is an N by N upper triangular binary matrix of the calculated logic
% domain
% MODES is an N by 1 column vector of completion modes 
% SST is an N by 1 column vector of scheduled start time 
%----------------
%Input:
%chromosome: row vector of decision values
%----------------
%Usage:
%PSM=updatepemdr(chromosome)
%----------------
%Example:
%This is is not executed as a stand-alone fuinction. Instead of specifying
%global values this function is cannot be run.
%----------------
%Prepositions and Requirements:
%1.)Global variable P must be specified.

function PSM=updatepemdr(chromosome)
global P %N by N upper triangular score matrix of logic domain
PEM=P;
nA=numel(P(P>0&P<1)); %Number of uncertainties and uncertain dependencies
n=numel(chromosome);  %Number of genes
N=(n-nA)/2;           %Number of tasks
if numel(PEM)>1
    g=1; %Gene g
    for C=1:size(PEM,1)
        for R=1:C
            if ((PEM(R,C)>0)&&(PEM(R,C)<1)) %Only uncertain values
                if numel(chromosome)>=g     %(=decision variable)
                    PEM(R,C)=chromosome(g); %is considered
                    g=g+1;
                end
            end
        end
    end
end
PSM=[round(PEM),round(chromosome(nA+1:nA+N))',chromosome(nA+N+1:end)'];