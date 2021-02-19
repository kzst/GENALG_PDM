%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Calculate an optimal project structure matrix and pareto-optimal resource
%allocation
%The algorithm is the implementation of Agile Project Management agent
%(APMa) by integer genetic algorithm and a multi-genetic algorithm to
%specify pareto-optimal resource allocation
%----------------
%Output:
%PSM: An n by N+2+nR+1 PSM=[DSM,TD,CD,RD,SST] matrix, where
% DSM is an N by N upper triangular binary matrix of the calculated logic
% domain
% TD is an N by 1 column vector of task durations 
% CD is an N by 1 column vector of cost demands
% RD is an N by nR matrix of resource demands
% SST is an N by 1 column vector of scheduled start time of tasks
%---------------- 
%Inputs:
%PDM: PDM=[PEM,TD,CD,RD] matrix (N by N+2+nR), 
%where PEM is an N by N upper triangular matrix of logic domain, 
% TD is an N by 1 column vector of task durations
% CD is an N by 1 column vector of cost demands
% RD is an N by nR matrix of resource demands
%const: 3 element vector of [Ct,Cc,Cs] values, where
% Ct is the time constraint (max constraint)
% Cc is the cost constraint (max constraint)
% Cs is the score constraint (min constraint)
%typefcn: Type of target function 
% 1=minTPT, 2=minTPC, 3=maxTPS, ~ composite
%---------------- 
%Usage:
%PSM=apmgenpr(PDM,const,typefcn)
%---------------- 
%Example:
%PDM=[triu(rand(10)*0.5+0.5),20*rand(10,1),30*rand(10,1),5*rand(10,3)];
%const=[percentt(PDM,1,.9),percentc(PDM,1,.9),percents(PDM,.7)];
%typefcn=999; %let target function be a composite target function
%tic;PSM=apmgenpr(PDM,const,typefcn);toc
%---------------- 
%Prepositions and Requirements:
%1.)The logic domian must be an upper triangular matrix, where the matrix
%elements are between 0 to 1 interval.
%2.)At least one matrix element is lower than 1 and upper than 0. => There
%are at least one decision variable
%3.)The elements of Time/Cost/Resource Domains are positive real numbers.

function PSM=apmgenpr(PDM,const,typefcn)
global PEM P Q T C R Ct Cc Cs Typ
%global values: 
% PEM is the input logic domain
% P is the score matrix of inclusion task dependencies/completions 
%the logic domain (initially PEM=P)
% Q is the score matrix of exclusion task dependencies/completions. 
%In this simulation Q=1-P
% T is the column vector of time demands (=the time domain)
% C is the column vector of cost demands (=the cost domain)
% R is the matrix of resource demands (=the resource domain)
% Ct is the time constraint
% Cc is the cost constraint
% Cs is the score constraint
% Typ is the selection number of the target functions (see above)

%---Initialization of the global values---

Typ=typefcn;
Ct=const(1);
Cc=const(2);
Cs=const(end);
N=size(PDM,1); %The number of activities 
pem=triu(PDM(:,1:N)); %Only the upper triangle is considered
T=PDM(:,N+1); %N+1-th column in PDM matrix is the time domain
C=PDM(:,N+2); %N+2-th column in PDM matrix is the cost domain
R=PDM(:,N+3:end); %N+3..N+3+nR-th columns in PDM matrix is the 
                  %resource domain
PEM=pem;
P=pem;
Q=ones(size(pem,1))-P; %Q=1-P

%---End of the initialization of the global values---

n=numel(PEM(PEM>0&PEM<1)); %number of uncertainties = 
% = number of variables in the genetic optimization
MinPopSize=100; %Minimal number of population
PopSize=min(max(1000,round((2^n)/100)),MinPopSize); %The population size 
%depends on the number of uncertainties

%---Set of genetic algorithm (GA)---

% PopulationType:bitstring = This is a binary genetic algorithm, where, the 
%every decision variables are 0 or 1.
% Display:off = There is no need any information to display
% PopulationSize:PopSize = The Population size in a parallel running is
%PopSize
% Vectorized:on and UseParallel:true parallelize the genetic algorithm
% Generations:100 and StallGenLimit:100 sets the maximal number of
%generations to 100
% TolCon:1e-6 means, that the GA will stop, if the given tolerance is
%reached
% CreationFcn:@createpem specifies createpem function when generating
%population.

options=gaoptimset('PopulationType','bitstring','Display','off',...
    'PopulationSize',PopSize,'Vectorized','on','UseParallel',true,...
    'Generations',100,'StallGenLimit',100,'TolCon',1e-6,...
    'CreationFcn', @createpem);

%---End of setting GA---

%---Runing GA---
% @targetfcn: the target function (=targetfcn)
% n: the number of variables
% options: The set of options (see above).
%Note: Because of the bitstring population, every genes are 0 or 1 => There 
%is no need upper and lower bounds. Nonlinear constrants are not allowed 
%when using bitstring population

chromosome=ga(@targetfcn,n,[],[],[],[],[],[],[],[],options);

%---End of runing GA---

PSM=updatepem(chromosome); %Now PSM only a bibary DSM matrix
PSM(diag(PSM)==0,:)=0; %All excluded tasks dependencies
PSM(:,diag(PSM)==0)=0; %is erased from the PSM matrix

%---Start of Pareto Optimizing---

SST=multires_solve(PSM,T,R); %Specify an SST for a Pareto-Optimal 
                    %resource allocation 

%---End of Pareto Optimizing---

%---Start of output formatting---

PSM=[PSM,T,C,R,SST];%The final PSM includes the original time/cost/resource 
                   %domains and the column vector of scheduled start time

%---End of output formatting---
                   
PSM(diag(PSM)==0,:)=0; %All excluded tasks dependencies/time/cost/resource
PSM(:,diag(PSM)==0)=0; %demands are erased from the PSM matrix


