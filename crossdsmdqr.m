%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Subfunction for Crossover of resource-constrained HDTQCTPs 
%----------------
%Output:
%chromosome=Recombined gene sequences
%---------------- 
%Inputs:
%DSM1: N by N binary upper triangle matrix of parent1's logic domain
%MODES1: N by 1 column vector of parent1's completion modes (=durations)
%SST1: N by 1 column vector of parent1's scheduled start time (SST)
%DSM2: N by N binary upper triangle matrix of parent2's logic domain
%MODES2: N by 1 column vector of parent2's completion modes (=durations)
%SST2: N by 1 column vector of parent2's scheduled start time (SST)
%---------------- 
%Usage:
%chromosome=crossdsmdqr(DSM1,MODES1,SST1,DSM2,MODES2,SST2)
%---------------- 
%Example:
%This is not executed as a stand-alone fuinction. Instead of specifying
%global values this function is cannot be run.
%---------------- 
%Prepositions and Requirements:
%1.)Global variable (P,T) must be specified.

function chromosome=crossdsmdqr(DSM1,MODES1,SST1,DSM2,MODES2,SST2)
global P T

%Global variables:
% P: N by N upper triangle score matrix of logic domain
% T: N by 1 column vector of task durations (time domain)
chromosome=[];
DSM1=round(DSM1); %DSM should contain binary values
DSM2=round(DSM2); %DSM should contain binary values
MODES1=round(MODES1); %MODES should contain discrete values
MODES2=round(MODES2); %MODES should contain discrete values
lout1=isfeasibledqr(DSM1,MODES1,SST1); %Feasibility check of parent1
lout2=isfeasibledqr(DSM2,MODES2,SST2); %Feasibility check of parent2
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
t=MODES;
for i=1:numel(MODES1)    
    if DSM(i,i)==0 %If task i is exluded, the completion mode is irrelevant
    else %The recombination of completion modes does not assume dominancy
        if rand<0.5 
            if DSM1(i,i)==0 %If only parent1's task i is exluded
                if DSM2(i,i)==0 %=>the child's completion mode will be 
                else            %parent2's completion mode
                    MODES(i)=MODES2(i); 
                end
            else %If only parent2's task i is exluded =>the child's 
                MODES(i)=MODES1(i); %completion mode will be parent1's          
            end                     %completion mode
        else
            if DSM2(i,i)==0 %If only parent2's task i is exluded
                 if DSM1(i,i)==0 %=>the child's completion mode will be 
                 else            %parent1's completion mode
                     MODES(i)=MODES1(i); 
                 end
            else %If only parent1's task i is exluded =>the child's 
                MODES(i)=MODES2(i); %completion mode will be parent2's      
            end                     %completion mode
        end
        if MODES(i)>0
            t(i)=T(i,MODES(i));
        end
    end
    chromosome=[chromosome,MODES(i)];
end
alpha=rand;
t=MODES; %Calculation of proposed SST
[~,EST,~,LST]=tptfast(DSM,t); %SST should be between EST and LST
SST=EST;
for i=1:numel(t)
    if DSM(i,i)>0 %..if task i is included
        if ((SST1(i)>=EST(i))&&(SST1(i)<=LST(i))&&(SST2(i)...
                >=EST(i))&&(SST2(i)<=LST(i))) %... If SSTs are in [EST,LST]
            SST(i)=alpha*SST1(i)+(1-alpha)*SST2(i);
        else
            if ((SST1(i)>=EST(i))&&(SST1(i)<=LST(i)))
                SST(i)=SST1(i); %... If only parent1'SST is in [EST,LST]
            else
                if ((SST2(i)>=EST(i))&&(SST2(i)<=LST(i)))
                    SST(i)=SST2(i);%... If only parent2'SST is in [EST,LST]
                else %... If nobody'SST are in [EST,LST]
                    SST(i)=alpha*EST(i)+(1-alpha)*LST(i);
                end
            end
        end
    end
    chromosome=[chromosome,SST(i)];
end
