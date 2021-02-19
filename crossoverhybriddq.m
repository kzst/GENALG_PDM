%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Crossover function for HCTQCTPs and Pareto-Optimal Resource Allocation
%----------------
%Standard output:
%xoverKids=Set of the result of the recombination
%---------------- 
%Standard inputs:
%parents=Set of chromosomes of parents
%options=Set of GA options not used
%GenomeLength=The size of the chromosomes
%FitnessFcn,unused=Fitness function, unused = not used parameters
%thisPopulation=Set of the population
%---------------- 
%Usage:
%xoverKids = crossoverhybriddq(parents,options,GenomeLength,FitnessFcn,...
    %unused,thisPopulation)
%---------------- 
%Example:
%This is not executed as a stand-alone fuinction. Instead of specifying
%global values this function is cannot be run.

function xoverKids  = crossoverhybriddq(parents,options,GenomeLength,...
    FitnessFcn,unused,thisPopulation)
nKids = length(parents)/2;
% Allocate space for the kids
xoverKids = zeros(nKids,GenomeLength);

% To move through the parents twice as fast as thekids are
% being produced, a separate index for the parents is needed
index = 1;
% for each kid...
for i=1:nKids
    % get parents
    r1 = parents(index);
    index = index + 1;
    r2 = parents(index);
    index = index + 1;
    % Decoding parents' chromosomes, which contains results of the
    % deciosions of uncertain tasks/dependencies (DSM) and completion modes
    % (MODES)
    PSM1=updatepemd(thisPopulation(r1,:));
    DSM1=PSM1(:,1:size(PSM1,1));
    MODES1=PSM1(:,end);
    PSM2=updatepemd(thisPopulation(r2,:));
    DSM2=PSM2(:,1:size(PSM2,1));
    MODES2=PSM2(:,end);
    % Randomly select half of the genes from each parent
    % This loop may seem like brute force, but it is twice as fast as the
    % vectorized version, because it does no allocation.
    xoverKids(i,:)=crossdsmdq(DSM1,MODES1,DSM2,MODES2);
end
