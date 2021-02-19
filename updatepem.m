%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Subfunction for conversion between chromosome and DSM
%----------------
%Output:
%DSM: N by N upper triangular binary matrix of logic domain
%----------------
%Input:
%chromosome: binary row vector of decision values
%----------------
%Usage:
%DSM=updatepem(chromosome)
%----------------
%Example:
%This is is not executed as a stand-alone fuinction. Instead of specifying
%global values this function is cannot be run.
%----------------
%Prepositions and Requirements:
%1.)Global variable P must be specified.

function DSM=updatepem(chromosome)
global P %N by N upper triangular score matrix of logic domain
PEM=P;
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
DSM=round(PEM); %DSM must be a binary upper triangular matrix!