%Author: Zsolt T. Koszty�n Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Subfunction for Crossover of HCTQCTPs and Pareto-Optimal Resource 
%Allocation
%----------------
%Output:
%chromosome=Recombined gene sequences
%---------------- 
%Inputs:
%DSM1: N by N binary upper triangle matrix of parent1's logic domain
%MODES1: N by 1 column vector of parent1's completion modes (=durations)
%DSM2: N by N binary upper triangle matrix of parent2's logic domain
%MODES2: N by 1 column vector of parent2's completion modes (=durations)
%---------------- 
%Usage:
%chromosome=crossdsmcq(DSM1,MODES1,DSM2,MODES2)
%---------------- 
%Example:
%This is not executed as a stand-alone fuinction. Instead of specifying
%global values this function is cannot be run.
%---------------- 
%Prepositions and Requirements:
%1.)Global variable P must be specified.

function chromosome=crossdsmcq(DSM1,MODES1,DSM2,MODES2)
global P % N by N upper triangle score matrix of logic domain
chromosome=[];
lout1=isfeasiblecq(DSM1,MODES1); %Feasibility check of parent1
lout2=isfeasiblecq(DSM2,MODES2); %Feasibility check of parent2
crit=.5; %The dominance factor. If crit=0.5, there is no dominant parent
if (lout1==lout2)
else
    if lout1==true %If parent1 is feasible, but parent2 not=>
        crit=.9;   %parent1 will be dominant
    else
        crit=.1;   %If parent2 is feasible, but parent1 not=>
    end            %parent2 will be dominant
end
for C=1:size(DSM1,1) %Calculating the logic domain for the child
    for R=1:C
        if ((P(R,C)>0)&&(P(R,C)<1)) 
            if rand<crit %Dominant genes are more probable inhereted
                if ((DSM1(R,R)>0)&&(DSM1(C,C)>0))
                    chromosome=[chromosome,DSM1(R,C)];
                else
                    chromosome=[chromosome,0]; %
                end
            else
                if ((DSM2(R,R)>0)&&(DSM2(C,C)>0))
                    chromosome=[chromosome,DSM2(R,C)];
                else
                    chromosome=[chromosome,0];
                end
            end        
        end
    end
end
DSM=updatepem(chromosome); %DSM is calculated, according to the child's 
                           %chromosomes
MODES=zeros(size(DSM,1),1);
alpha=rand;
for i=1:numel(MODES1)
    if DSM(i,i)==0 %If task i is exluded, the completion mode is irrelevant
    else
        if MODES1(i)==0 %If parent1's task i is exluded, but parent2's not
            if MODES2(i)==0 %=>the child's completion mode will be 
            else %parent2's completion
                MODES(i)=MODES2(i);
            end
        else        %If parent2's task i is exluded, but parent1's not
            if MODES2(i)==0 %=>the child's completion mode will be 
                MODES(i)=MODES1(i); %parent1's completion
            else    %If both parents tasks are included, the child's mode
                MODES(i)=alpha*MODES1(i)+(1-alpha)*MODES2(i); %are randomly
            end                                                   %selected
        end
    end
    chromosome=[chromosome,MODES(i)]; 
end