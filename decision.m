%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Subfunction for conversion between chromosome and DSM
%----------------
%Output:
%chromosome=row vector of decision values
%---------------- 
%Input:
%PEM: N by N upper triangular matrix 
%---------------- 
%Usage:
%chromosome=decision(PEM)
%---------------- 
%Example:
%This is not executed as a stand-alone fuinction. Instead of specifying
%global values this function is cannot be run.
%---------------- 
%Prepositions and Requirements:
%1.)Global variable P must be specified.

function chromosome=decision(PEM)
global P %N by N upper triangular score matrix of logic domain
chromosome=[];
PEM=triu(PEM); %Only the upper triangular (sub)matrix is considered
if numel(PEM)==1
    if PEM==P
    else
        chromosome=PEM;
    end
else
    if numel(PEM)>1
        for C=1:size(PEM,1)
            for R=1:C
                if ((P(R,C)>0)&&(P(R,C)<1)) %Only uncertain values 
                   chromosome=[chromosome,PEM(R,C)]; %(=decision variable) 
                end                                  %is considered
            end
        end
    end
end