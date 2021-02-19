%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Calculate an optimal project structure matrix and optimal resource
%allocation for a resource-constraint hybrid discrete time-quality-cost
%trade-off problem (RCHDTQCTP)
%The algorithm is the implementation of resource-constrained Hybrid Project 
%Management agent considering quality parameters (rcHPMa(q)) by hybrid 
%genetic algorithm 
%----------------
%Output:
%PSM: An N by N+3+nR+1 PSM=[DSM,TD,CD,QD,RD,SST] matrix, where
% DSM is an N by N upper triangular binary matrix of the calculated logic
% domain
% TD is an N by 1 column vector of task durations 
% CD is an N by 1 column vector of cost demands
% QD is an N by 1 column vector of quality demands
% RD is an N by nR matrix of resource demands
% SST is an N by 1 column vector of scheduled start time of tasks
%---------------- 
%Inputs:
%PDM: PDM=[PEM,TD,CD,QD,RD] matrix (N by N+(3+nR)w), 
%where PEM is an N by N upper triangular matrix of logic domain, 
% TD is an N by w matrix of task durations
% CD is an N by w matrix of cost demands
% QD is an N by w matrix of quality parameters
% RD is an N by w*nR matrix of resource demands
%const: 3+nR+1 elements row vector of [Ct,Cc,Cq,CR,Cs] values, where
% Ct is the time constraint (max constraint)
% Cc is the cost constraint (max constraint)
% Cq is the quality constraint (min constraint)
% CR are the resource constraints (max constraint)
% Cs is the score constraint (min constraint)
%typetfcn: Type of target function 
% 0=maxTPQ, 1=minTPT, 2=minTPC, 3=maxTPS, 4=minUF, ~ composite
%---------------- 
%Usage:
%PSM=hpmgendqr(PDM,const,typefcn)
%---------------- 
%Example:
%PDM=[triu(rand(10)*.5+.5),20*rand(10,3),30*rand(10,3),rand(10,3),...
%   5*rand(10,6)];
%const=[percentt(PDM,3,.9),percentc(PDM,3,.9),percentqr(PDM,3,1)',...
%   percentq(PDM,3,.7),percents(PDM,.7)];
%typefcn=999; %let target function be a composite target function
%tic;PSM=hpmgendqr(PDM,const,typefcn);toc
%---------------- 
%Prepositions and Requirements:
%1.)The logic domian must be an upper triangular matrix, where the matrix
%elements are between 0 to 1 interval.
%2.)At least one matrix element is lower than 1 and upper than 0. => There
%are at least one decision variable
%3.)The number of modes(w) is a positive integer
%4.)The number of resources (nR) is a positive number
%5.)The elements of Time/Cost/Quality/Resource Domains are positive real 
%numbers.
%6.)Usually a monotonity are assumed, which, means: if tk,i<tk,j => 
%ck,i>=ck,j and qk,i<=qk,j, where 1<=k<=N is the k-th task, 1<=i,j<=w are 
%the selected modes. Nevertheless, this assumption is not required in this 
%simulation. 

function PSM=hpmgendqr(PDM,const,typefcn)
global PEM P Q T C q R Ct Cc CR Cq Cs Typ QD

%global values: 
% PEM is the input logic domain
% P is the score matrix of inclusion task dependencies/completions 
%the logic domain (initially PEM=P)
% Q is the score matrix of exclusion task dependencies/completions. 
%In this simulation Q=1-P
% T is the N by w matrix of time demands (=the time domain)
% C is the N by w matrix of of cost demands (=the cost domain)
% q is the N by w matrix of of quality parameters (=the quality domain)
% R is the N by w*nR matrix of of resource demands (=the resource domain)
% Ct is the time constraint
% Cc is the cost constraint
% Cq is the quality constraint
% CR is the row vector resource constraint
% Cs is the score constraint
% Typ is the selection number of the target functions (see above)

%---Initialization of the global values---

Typ=typefcn;
Ct=const(1);
Cc=const(2);
Cq=const(3);
CR=const(4:end-1);
Cs=const(end);

N=size(PDM,1); %The number of activities 
M=size(PDM,2); %M=N+(3+nR)w
nR=numel(const)-4; %Number of resources
w=(M-N)/(nR+3); %Number of modes

pem=triu(PDM(:,1:N)); %only the upper triangle will be considered

T=PDM(:,N+1:N+w); %N+1..N+w-th columns in PDM matrix is the time domain
C=PDM(:,N+1+w:N+2*w); %N+w+1..N+2*w-th columns in PDM matrix is the cost 
                      %domain
q=PDM(:,N+1+2*w:N+3*w); %N+2*w+1..N+3*w-th columns in PDM matrix is the
                        %quality domain
QD=q;
R=PDM(:,N+1+3*w:end); %N+3*w+1..N+(3+nR)*w-th columns in PDM matrix is the
                      %resource domain
                      
PEM=pem;
P=pem;
Q=ones(N)-P; %Q=1-P

%---End of the initialization of the global values---

n=numel(PEM(PEM>0&PEM<1)); %Number of uncertainties

MinPopSize=100; %Minimal number of population
PopSize=min(max(1000,round((2^n)/100)),MinPopSize); %The population size 
%depends on the number of uncertainties

%---Set of genetic algorithm (GA)---

chromosomedqr=initpopdqr(PopSize); %Specify the initial population

% Display:off = There is no need any information to display 
% PopulationSize:PopSize = The Population size in a parallel running is 
%PopSize 
% Vectorized:on and UseParallel:true parallelize the genetic algorithm 
% Generations:100 and StallGenLimit:100 sets the maximal number of 
%generations to 100 
% TolCon:1e-6 means, that the GA will stop, if the given tolerance is 
%reached 
% InitialPopulation:chromosomedqr chromosomedqr predefines set of chromosomes 
%are used as an initial population
% CreationFcn:@createpemdqr specifies createpemdqr function when generating 
%population. 
% HybridFcn:@fmincon specifies a hybrid function which is used after 
%running GA. 
% CrossoverFcn:@crossoverhybriddqr specifies a unique crossoverhybriddqr  
%crossover function, which recombines two PSM=[DSM,MODES,SST] matrices   
% MutationFcn:@mutationadaptfeasible randomly generates directions that are 
%adaptive with respect to the last successful or unsuccessful generation 
% SelectionFcn:{@selectiontournament,4} Tournament selection chooses each  
%parent by choosing Tournament size players at random and then choosing the 
%best individual out of that set to be a parent. Tournament size is 4. 

options=gaoptimset('Display','off','PopulationSize',PopSize,...
    'Vectorized','on','UseParallel',true,...
    'Generations',100,'StallGenLimit',100,'TolCon',1e-6,...
    'InitialPopulation',chromosomedqr,'CreationFcn',@createpemdqr,...
    'CrossoverFcn',@crossoverhybriddqr,...
    'MutationFcn',@mutationadaptfeasible,...
    'SelectionFcn',{@selectiontournament,4});

%---Setting the global lower/upper bounds---

[~,EST]=tptfast(floor(PEM),min(T,[],2)); %EST when all uncertain 
                                  %realizations are excluded from the DSM
[~,~,~,LST]=tptfast(ceil(PEM),max(T,[],2)); %LST when all uncertain
                                  %realizations are included from the DSM
LST(LST<0)=0;
LB=[zeros(1,numel(PEM(PEM>0&PEM<1))),ones(1,N),EST'];
UB=[ones(1,numel(PEM(PEM>0&PEM<1))),ones(1,N)*w,LST'];

%---End of bounds settings--- 

%---End of setting GA---  

%---Runing GA--- 

% @targetfcndqr: the target function (=targetfcndqr) 
% n+2*N: the number of variables (n=number of uncertainties + 2*N=2*number  
%of activities. The first n variable is a binary decision variable. The  
%second part of chromosome(n+1..n+N) are integer values and contains modes,  
%the last part of chromosome(n+N+1..n+2*N) contains the scheduled start 
%time of activities.
% LB,UB: Lower/Upper Bounds 
% @constrainrtpardqr: The nonlinear constraint function 
% options: The set of options (see above). 
%Note: Despite binary,integer and real values have to be seeked, in order  
%to apply custom creation, crossover and hybrid functions the binary and  
%integer values are considered as rounded continuous values 

chromosome=ga(@targetfcndqr,n+2*N,[],[],[],[],LB,UB,@constrainrtpardqr,options);

%---End of runing GA---

%---Start of output formatting---   

PSM=updatepemdr(chromosome);%PSM=[DSM,T,SST] = Logic Domain and the vector 
                            %of the selected modes and the vector of
                            %scheduled start time of activities

PSM(diag(PSM)==0,:)=0; %All excluded tasks dependencies/time demands/  
PSM(:,diag(PSM)==0)=0; %cost demands are erased from the PSM matrix                               
                                
PSM=pmtopsmqr(PSM,T,C,q,R); %The output PSM = [DSM,TD,CD,QD,RD,SST] 

PSM(diag(PSM)==0,:)=0; %All excluded tasks dependencies/time demands/   
PSM(:,diag(PSM)==0)=0; %cost demands are erased from the PSM matrix         

%---End of output formatting---