%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Calculate an optimal project structure matrix and optimal resource
%allocation for a hybrid (continouos) multi-objective multi-mode resource 
%constraint project scheduling problem with quality ((c)HMOMRCPSP(q)) by 
%hybrid genetic algorithm 
%----------------
%Output:
%PSM: An N by N+2+nR+1 PSM=[DSM,TD,CD,RD,SST] matrix, where
% DSM is an N by N upper triangular binary matrix of the calculated logic
% domain
% TD is an N by 1 column vector of tasknd optimal resource
%allocation for resource constraint hybrid time-cost trade-off problem
%(RCHDTCTP) durations 
% CD is an N by 1 column vector of cost demands
% RD is an N by nR matrix of resource demands
% SST is an N by 1 column vector of scheduled start time of tasks
%---------------- 
%Inputs:
%PDM: PDM=[PEM,TD,CD,RD] matrix (N by N+(2+nR)w), 
%where PEM is an N by N upper triangular matrix of logic domain, 
% TD is an N by w matrix of task durations
% CD is an N by w matrix of cost demands
% RD is an N by w*nR matrix of resource demands
%const: 2+nR+1 elements row vector of [Ct,Cc,CR,Cs] values, where
% Ct is the time constraint (max constraint)
% Cc is the cost constraint (max constraint)
% CR are the resource constraints (max constraint)
% Cs is the score constraint (min constraint)
%---------------- 
%Usage:
%PSM=parhpmgendr(PDM,const)
%---------------- 
%Example:
%PDM=[triu(rand(10)*.5+.5),20*rand(10,3),30*rand(10,3),5*rand(10,6)];
%const=[percentt(PDM,3,.9),percentc(PDM,3,.9),percentr(PDM,3,1)',...
%  percents(PDM,.7)];
%typefcn=999; %let target function be a composite target function
%tic;PSM=parhpmgendr(PDM,const);toc
%---------------- 
%Prepositions and Requirements:
%1.)The logic domian must be an upper triangular matrix, where the matrix
%elements are between 0 to 1 interval.
%2.)At least one matrix element is lower than 1 and upper than 0. => There
%are at least one decision variable
%3.)The number of modes(w) is a positive integer
%4.)The number of resources (nR) is a positive number
%5.)The elements of Time/Cost/Resource Domains are positive real numbers.
%6.)Usually a monotonity are assumed, which, means: if tk,i<tk,j => 
%ck,i>=ck,j, where 1<=k<=N is the k-th task, 1<=i,j<=w are the selected
%modes. Nevertheless, this assumption is not required in this simulation. 

function PSM=parhpmgendr(PDM,const)
global PEM P Q T C R Ct Cc CR Cs 

%global values: 
% PEM is the input logic domain
% P is the score matrix of inclusion task dependencies/completions 
%the logic domain (initially PEM=P)
% Q is the score matrix of exclusion task dependencies/completions. 
%In this simulation Q=1-P
% T is the N by w matrix of time demands (=the time domain)
% C is the N by w matrix of cost demands (=the cost domain)
% R is the N by w*nR matrix of resource demands (=the resource domain)
% Ct is the time constraint
% Cc is the cost constraint
% CR are 1 by nR row vector of resource constraints
% Cs is the score constraint

%---Initialization of the global values---

Ct=const(1);
Cc=const(2);
CR=const(3:end-1);
Cs=const(end);
N=size(PDM,1); %The number of activities 
M=size(PDM,2); %M=N+(2+nR)w
nR=numel(const)-3; %Number of resources
w=(M-N)/(nR+2); %Number of nodes
pem=triu(PDM(:,1:N)); %only the upper triangle will be considered

T=PDM(:,N+1:N+w); %N+1..N+w-th columns in PDM matrix is the time domain
C=PDM(:,N+1+w:N+2*w); %N+w+1..N+2*w-th columns in PDM matrix is the cost 
                      %domain
R=PDM(:,N+1+2*w:end); %N+2*w+1..N+(2+nR)*w-th columns in PDM matrix is the
                      %resource domain
PEM=pem;
P=pem;
Q=ones(N)-P; %Q=1-P

%---End of the initialization of the global values---

n=numel(PEM(PEM>0&PEM<1));

MinPopSize=100; %Minimal number of population
PopSize=min(max(1000,round((2^n)/100)),MinPopSize); %The population size 
%depends on the number of uncertainties

%---Set of genetic algorithm (GA)---

chromosomedr=initpopdr(PopSize); %Specify the initial population

% Display:off = There is no need any information to display
% PopulationSize:PopSize = The Population size in a parallel running is
%PopSize
% Vectorized:on and UseParallel:true parallelize the genetic algorithm
% Generations:100 and StallGenLimit:100 sets the maximal number of
%generations to 100
% TolCon:1e-6 means, that the GA will stop, if the given tolerance is
%reached
% InitialPopulation:chromosomedr chromosomedr predefines set of chromosomes 
%are used as an initial population
% CreationFcn:@createpemdr specifies createpemdr function when generating
%population.
% HybridFcn:@fmincon specifies a hybrid function which is used after
%running GA.
% CrossoverFcn:@crossoverhybriddr specifies a unique crossoverhybriddr 
%crossover function, which recombines two PSM=[DSM,MODES,SST] matrices  
% MutationFcn:@mutationadaptfeasible randomly generates directions that are
%adaptive with respect to the last successful or unsuccessful generation
% SelectionFcn:{@selectiontournament,4} Tournament selection chooses each 
%parent by choosingï¿½Tournament sizeï¿½players at random and then choosing the
%best individual out of that set to be a parent.ï¿½Tournament sizeï¿½isï¿½4.

options=gaoptimset('Display','off','PopulationSize',PopSize,...
    'Vectorized','on','UseParallel',true,...
    'Generations',100,'StallGenLimit',100,'TolCon',1e-6,...
    'InitialPopulation',chromosomedr,'CreationFcn',@createpemdr,...
    'CrossoverFcn',@crossoverhybriddr,...
    'MutationFcn',@mutationadaptfeasible,...
    'SelectionFcn',{@selectiontournament,4});

%---Setting the global lower/upper bounds---

[~,EST]=tptfast(floor(PEM),min(T,[],2)); %EST when all uncertain 
                                  %realizations are excluded from the DSM
[~,~,~,LST]=tptfast(ceil(PEM),max(T,[],2)); %LST when all uncertain
                                  %realizations are included from the DSM
LST(LST<0)=0;
LB=[zeros(1,numel(PEM(PEM>0&PEM<1))),ones(1,N),EST'];          %Lower Bound
UB=[ones(1,numel(PEM(PEM>0&PEM<1))),ones(1,N)*w,LST'];         %Upper Bound
for i=1:numel(LB)
    if LB(i)>UB(i)
        UB(i)=LB(i);
    end
end
%---End of bounds settings---

%---End of setting GA---

%---Runing GA---
% @targetfcndr: the target function (=targetfcndr)
% n+2*N: the number of variables (n=number of uncertainties + 2*N=2*number 
%of activities. The first n variable is a binary decision variable. The 
%second part of chromosome(n+1..n+N) are integer values and contains modes, 
%the last part of chromosome(n+N+1..n+2*N) contains the scheduled start
%time of activities.
% LB,UB: Lower/Upper Bounds
% @constrainrtpardr: The nonlinear constraint function
% options: The set of options (see above).
%Note: Despite binary,integer and real values have to be seeked, in order 
%to apply custom creation, crossover and hybrid functions the binary and 
%integer values are considered as rounded continuous values

X=gamultiobj(@targetfcndrpar,n+2*N,[],[],[],[],LB,UB,@constrainrtpardr,...
    options);

chromosome=X(1,:);
%---End of runing GA---  
 
%---Start of output formatting---  

PSM=updatepemdr(chromosome);%PSM=[DSM,T,SST] = Logic Domain and the vector 
                            %of the selected modes and the vector of
                            %scheduled start time of activities
                            
PSM(diag(PSM)==0,:)=0; %All excluded tasks dependencies/time demands/ 
PSM(:,diag(PSM)==0)=0; %cost demands are erased from the PSM matrix

PSM=pmtopsmr(PSM,T,C,R); %The output PSM = [DSM,TD,CD,RD,SST]

PSM(diag(PSM)==0,:)=0; %All excluded tasks dependencies/time demands/ 
PSM(:,diag(PSM)==0)=0; %cost demands are erased from the PSM matrix

%---End of output formatting---