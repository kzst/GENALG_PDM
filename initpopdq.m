%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Create initial population for HDTCTPs and HDTQCTPs and Pareto-Optimal 
%Resource Allocations
%----------------
%Output:
%chromosome=PopSize by n+N matrix of included/excluded dependencies/tasks
%and completion modes
%---------------- 
%Input:
%PopSize=Size of the population
%---------------- 
%Usage:
%chromosome=initpopdq(PopSize)
%---------------- 
%Example:
%This is not executed as a stand-alone fuinction. Instead of specifying
%global values this function is cannot be run.
%---------------- 
%Prepositions and Requirements:
%1.)Global values P,T must be specified.

function chromosome=initpopdq(PopSize)
global P T
%global values: 
% P is the score matrix of inclusion task dependencies/completions 
% T is the column vector of time demands (=the time domain)
nA=sum(diag(P)>0&diag(P)<1); %Number of uncertain acitivieis
n=numel(P(P>0&P<1)); %Number of uncertain activities+uncertain dependencies
w=size(T,2); %Number of modes
N=size(P,1); %Number of activities
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
