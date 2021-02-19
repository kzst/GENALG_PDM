%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Calculate an optimal project structure matrix
%The algorithm is the implementation of Agile Project Management agent
%(APMa) by integer genetic algorithm
%----------------
%Output:
%PSM: An n by N+4 PSM=[DSM,TD,CD,QD,EST] matrix, where
% DSM is an N by N upper triangular binary matrix of the calculated logic
% domain
% TD is an N by 1 column vector of task durations 
% CD is an N by 1 column vector of cost demands
% QD is an N by 1 column vector of quality parameters
% EST is an N by 1 column vector of early start time of tasks
%---------------- 
%Inputs:
%PDM: PDM=[PEM,TD,CD,QD] matrix (n by n+3), 
%where PEM is an N by N upper triangular matrix of logic domain, 
% TD is an N by 1 column vector of task durations
% CD is an N by 1 column vector of cost demands
% QD is an N by 1 column vector of quality parameters
%const: 4 element vector of [Ct,Cc,Cq,Cs] values, where
% Ct is the time constraint (max constraint)
% Cc is the cost constraint (max constraint)
% Cq is the quality constraint (min constraint)
% Cs is the score constraint (min constraint)
%typefcn: Type of target function 
% 0=maxTPQ, 1=minTPT, 2=minTPC, 3=maxTPS, ~ composite
%---------------- 
%Usage:
%PSM=apmgenq(PDM,const,typefcn)
%---------------- 
%Example:
%PDM=[triu(rand(10)*.5+.5),20*rand(10,1),30*rand(10,1),rand(10,1)];
%const=[percentt(PDM,1,.9),percentc(PDM,1,.9),percentq(PDM,1,.7),percents(PDM,.7)];
%typefcn=999; %let target function be a composite target function
%tic;PSM=apmgenq(PDM,const,typefcn);toc
%---------------- 
%Prepositions and Requirements:
%1.)The logic domian must be an upper triangular matrix, where the matrix
%elements are between 0 to 1 interval.
%2.)At least one matrix element is lower than 1 and upper than 0. => There
%are at least one decision variable
%3.)The elements of Time/Cost/Quality Domains are positive real numbers.

function PSM=apmgenq(PDM,const,typefcn)
global PEM P Q T C q Ct Cc Cq Cs Typ q
%global values: 
% PEM is the input logic domain
% P is the score matrix of inclusion task dependencies/completions 
%the logic domain (initially PEM=P)
% Q is the score matrix of exclusion task dependencies/completions. 
%In this simulation Q=1-P
% T is the column vector of time demands (=the time domain)
% C is the column vector of cost demands (=the cost domain)
% q is the column vector of quality parameters (=the quality domain)
% Ct is the time constraint
% Cc is the cost constraint
% Cq is the quality constraint
% Cs is the score constraint
% Typ is the selection number of the target functions (see above)

%---Initialization of the global values---

Typ=typefcn;
Ct=const(1);
Cc=const(2);
Cq=const(3);
Cs=const(end);
N=size(PDM,1); %The number of activities 
pem=triu(PDM(:,1:N)); %Only the upper triangle is considered
T=PDM(:,N+1); %N+1-th column in PDM matrix is the time domain
C=PDM(:,N+2); %N+2-th column in PDM matrix is the cost domain
q=PDM(:,N+3); %N+3-th column in PDM matrix is the quality domain
QD=q;
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
% @targetfcnq: the target function (=targetfcnq)
% n: the number of variables
% options: The set of options (see above).
%Note: Because of the bitstring population, every genes are 0 or 1 => There 
%is no need upper and lower bounds. Nonlinear constrants are not allowed 
%when using bitstring population

chromosome=ga(@targetfcnq,n,[],[],[],[],[],[],[],[],options);

%---End of runing GA---

%---Start of output formatting---

PSM=updatepem(chromosome); %Now PSM only a bibary DSM matrix
[~,EST]=tptfast(PSM,T); %Determine EST
PSM=[PSM,T,C,q,EST];%The final PSM includes the original time/cost/quaility 
                   %domains and the column vector of earlies start time
                   
%---End of output formatting---

PSM(diag(PSM)==0,:)=0; %All excluded tasks dependencies/time/cost/quality
PSM(:,diag(PSM)==0)=0; %demands are erased from the PSM matrix