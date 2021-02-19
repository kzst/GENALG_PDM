%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Create population for Agile Project Scheduling.Problems
%----------------
%Standard output:
%Population=PopSize by nA binary matrix of included/excluded
%dependencies/tasks
%---------------- 
%Standard inputs:
%GenomeLength=nA number of uncertainties (not used)
%FittnessFcn The fitness function (not used)
%options options of GA (not used)
%---------------- 
%Usage:
%Population=createpem(GenomeLength, FitnessFcn, options)
%---------------- 
%Example:
%This is not executed as a stand-alone fuinction. Instead of specifying
%global values this function is cannot be run.
%---------------- 
%Prepositions and Requirements:
%1.)Global variable P must be specified.

function Population=createpem(GenomeLength, FitnessFcn, options)
global P %The N by N score matrix of task/dependency inclusion
nA=sum(diag(P)>0&diag(P)<1); %The number of uncertainties
PopSize=options.PopulationSize; %The number of population
chromosomediag=round(rand(PopSize,nA)); 
chromosome=[];
for i=1:PopSize
    chromosome=[chromosome;randoutdiag(P,chromosomediag(i,:))];
    DSM=updatepem(chromosome(i,:));
    DSM(diag(DSM)==0,:)=0; %Exclude dependencies of excluded tasks
    DSM(:,diag(DSM)==0)=0; %Exclude dependencies of excluded tasks
    chromosome(i,:)=decision(DSM);%chromsomes contains uncertain
end          %realizations
Population=chromosome;