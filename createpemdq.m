%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Create population for HDTCTPs and HDTQCTPs and Pareto-Optimal 
%Resource Allocations
%----------------
%Standard output:
%Population=PopSize by n+N matrix of included/excluded dependencies/tasks
%and completion modes
%---------------- 
%Standard inputs:
%GenomeLength=n+N number of uncertainties + number of tasks (not used)
%FittnessFcn The fitness function (not used)
%options options of GA (not used)
%---------------- 
%Usage:
%Population=createpemdq(GenomeLength, FitnessFcn, options)
%---------------- 
%Example:
%This is not executed as a stand-alone fuinction. Instead of specifying
%global values this function is cannot be run.
%---------------- 
%Prepositions and Requirements:
%1.)Global values P,T must be specified.

function Population= createpemdq(GenomeLength, FitnessFcn, options)
global P T
%global values: 
% P is the score matrix of inclusion task dependencies/completions 
% T is the column vector of time demands (=the time domain)
nA=sum(diag(P)>0&diag(P)<1); %Number of uncertain acitivieis
n=numel(P(P>0&P<1)); %Number of uncertain activities+uncertain dependencies
w=size(T,2); %Number of modes
N=size(P,1); %Number of activities
PopSize=options.PopulationSize; 
chromosomediag=round(rand(PopSize,nA));
chromosome=zeros(PopSize,n+N);
for I=1:PopSize
    uncertainties=randoutdiag(P,chromosomediag(I,:)); %First, the set of 
    DSM=updatepem(uncertainties);    %completed activities are specified
    DSM(diag(DSM)==0,:)=0; %Exclude dependencies of excluded tasks
    DSM(:,diag(DSM)==0)=0; %Exclude dependencies of excluded tasks
    MODES=ceil(rand(1,N)*w); %...between in {1,w}
    MODES(diag(DSM)==0)=0; %A duration is set to be zero when a task is 
                           %excluded
    chromosome(I,:)=[decision(DSM),MODES]; %chromsomes contains uncertain
end          %realizations and completion modes = task durations
Population=chromosome;