%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Calculate the constraints for the hybrid discrete time-cost trade-off 
%problems
%----------------
%Outputs:
%c is a vector of inequities.
%ceq is null vector of Eq.s
%---------------- 
%Input:
%chromosome: set of populations 
%---------------- 
%Usage:
%[c,ceq]=constrainrtpard(chromosome)
%---------------- 
%Example:
%This is not executed as a stand-alone fuinction. Instead of specifying
%global values this function is cannot be run.
%---------------- 
%Prepositions and Requirements:
%1.)Global values P,Q,T,C,Ct,Cc,Cs must be specified.
%2.)P should be an upper triangular matrix of the logic domain.

function [c,ceq]=constrainrtpard(chromosome)
global P Q T C Ct Cc Cs
%global values: 
% P is the PopSize by PopSize scores matrix of task/dependency inclusion
% Q is the PopSize by PopSize scores matrix of task/dependency exclusion
% T is the column vector of time demands (=the time domain)
% C is the column vector of cost demands (=the cost domain)
% Ct is the time constraint
% Cc is the cost constraint
% Cs is the score constraint

chromosome=round(chromosome); %Every genes of the chromosome must be an
                              %integer value
PopSize=numel(chromosome(:,1)); %The size of the population
N=size(P,1); %Number of activities
c=zeros(PopSize,3); %Initialize the PopSize by 3 matrix of inequities
tpts=zeros(PopSize,1); %Initialize the vector of TPTs
tpcs=zeros(PopSize,1); %Initialize the vector of TPCs

for I=1:PopSize
    PSM=updatepemd(chromosome(I,:));%PSM=[DSM,MODES] = Logic Domain and the 
                                     %vector of the selected modes
    DSM=PSM(:,1:N); %DSM is a binary upper triangular matrix
    MODES=PSM(:,end); %MODES is a vector of selected modes
    TD=zeros(N,1);
    CD=zeros(N,1);
    for i=1:N
        if MODES(i)==0
            TD(i)=0;
            CD(i)=0;
        else
            TD(i)=T(i,MODES(i));
            CD(i)=C(i,MODES(i));
        end
    end
    tpts(I)=tptfast(DSM,TD); %Calculate TPT for the member I of the 
                             %population
    tpcs(I)=diag(DSM)'*CD;   %Calculate TPC for the member I of the 
                             %population
    c(I,1)=tpts(I)-Ct; %TPT should be lower than the time constraint
    c(I,2)=tpcs(I)-Cc; %TPC should be lower than the cost constraint
    c(I,3)=Cs-maxscore_PEM(DSM,P,Q); %TPS should be greater than Cs
end
ceq=[];