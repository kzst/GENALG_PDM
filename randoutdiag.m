%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Generate random dependendencies (= Generate project structures)
%----------------
%Output:
%chromosome is a row vector of genes (=task
%inclusions/exclusions/dependencies)
%---------------- 
%Inputs:
%PEM is an N by N upper triangular matrix of logic domain
%diagchromosome is a row vector of decisions of task completions
%---------------- 
%Usage:
%chromosome=randoutdiag(PEM,diagchromosome)
%---------------- 
%Example:
%This is not executed as a stand-alone fuinction.

function chromosome=randoutdiag(PEM,diagchromosome)
global P
chromosome=[];
SNPM=updatepemdiag(diagchromosome);
if numel(PEM)==1
    if PEM==P
    else
        chromosome=PEM;
    end
else
    if numel(PEM)>1
        g=1;
        for C=1:size(PEM,1)
            for R=1:C
                if ((P(R,C)>0)&&(P(R,C)<1))
                    if R==C
                        if g<=numel(diagchromosome)
                            chromosome=[chromosome,diagchromosome(g)];
                            g=g+1;
                        end
                    else
                        if (SNPM(R,R)==1)&&(SNPM(C,C)==1) 
                            %Dependencies can occur if tasks, which has to
                            %be connected, decided as to be completed
                            chromosome=[chromosome,round(rand(1))];
                        else
                            chromosome=[chromosome,0];
                        end
                    end
                    
                end
            end
        end
    end
end