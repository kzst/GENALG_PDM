%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Calculate an optimal project structure matrix and pareto-optimal resource
%allocation for pareto-optimal resource balanced hybrid continouos
%time-quality-cost trade-off problem (PO-HCTQCTP)
%The algorithm is the implementation of Hybrid Project Management agent
%(HPMa) by hybrid genetic algorithm
%----------------
%Output:
%PSM: An N by N+3+nR+1 PSM=[DSM,TD,CD,QD,RD,SST] matrix, where
% DSM is an N by N upper triangular binary matrix of the calculated logic
% domain
% TD is an N by 1 column vector of task durations 
% CD is an N by 1 column vector of cost demands
% QD is an N by 1 column vector of quality parameters
% RD is an N by nR matrix of resource demands
% SST is an N by 1 column vector of scheduled start time of tasks
%---------------- 
%Inputs:
%PDM: PDM=[PEM,TD,CD,QD,RD] matrix (N by N+(3+nR)*2, 
%where PEM is an N by N upper triangular matrix of logic domain, 
% TD is an N by 2 matrix of task durations
% CD is an N by 2 matrix of cost demands
% QD is an N by 2 matrix of quality parameters
% RD is an N by 2*nR matrix of resource demands
%const: 3 element vector of [Ct,Cc,Cs] values, where
% Ct is the time constraint (max constraint)
% Cc is the cost constraint (max constraint)
% Cq is the quality constraint (min constraint)  
% Cs is the score constraint (min constraint)
%typetfcn: Type of target function 
% 0=maxTPQ, 1=minTPT, 2=minTPC, 3=maxTPS, ~ composite 
%---------------- 
%Usage:
%PSM=hpmgencqpr(PDM,const,typefcn)
%---------------- 
%Example:
%PDM=[triu(rand(10)*.5+.5),20*rand(10,2),30*rand(10,2),rand(10,2),...
%  5*rand(10,6)];
%const=[percentt(PDM,2,.9),percentc(PDM,2,.9),percentq(PDM,2,.7),...
%  percents(PDM,.7)];
%typefcn=999; %let target function be a composite target function
%tic;PSM=hpmgencqpr(PDM,const,typefcn);toc
%---------------- 
%Prepositions and Requirements:
%1.)The logic domian must be an upper triangular matrix, where the matrix
%elements are between 0 to 1 interval.
%2.)At least one matrix element is lower than 1 and upper than 0. => There
%are at least one decision variable
%3.)The number of modes(w) is 2,
%4.)The elements of Time and Cost Domains are positive real numbers.
%5.)Usually a monotonity are assumed, which, means: if tk,i<tk,j => 
%ck,i>=ck,j and qk,i<=qk,j, where 1<=k<=N is the k-th task, 1<=i,j<=w are 
%the selected modes. Therefore the algorithm match maximal cost/maximal 
%resource/minimal quality to the minimal durations and minimal cost/
%minimal resource/maximal quality to maximal durations. 

function PSM=hpmgencqpr(PDM,const,typefcn)
global PEM P Q T C q R Ct Cc Cq Cs Typ QD
%global values: 
% PEM is the input logic domain
% P is the score matrix of inclusion task dependencies/completions 
%the logic domain (initially PEM=P)
% Q is the score matrix of exclusion task dependencies/completions. 
%In this simulation Q=1-P
% T is the N by 2 matrix of time demands (=the time domain)
% C is the N by 2 matrix of cost demands (=the cost domain)
% q is the N by 2 matrix of quality demands (=the quality domain)  
% R is the N by 2*nR matrix of resource demands (=the resource domain)
% Ct is the time constraint
% Cc is the cost constraint
% Cq is the quality constraint
% Cs is the score constraint
% Typ is the selection number of the target functions (see above)

%---Initialization of the global values---

w=2; %In the case of CTCTP, the number of modes(w) is always 2. 
Typ=typefcn;
Ct=const(1);
Cc=const(2);
Cq=const(3);
Cs=const(end);

N=size(PDM,1); %The number of activities 
M=size(PDM,2); %M=N+3*2+2*nR
nR=(M-N-6)/2; %nR: the number of resources 
pem=triu(PDM(:,1:N)); %Only the upper triangle will be considered
TD=PDM(:,N+1:N+2);    %N+1..N+2-th columns in PDM matrix is the time domain
CD=PDM(:,N+3:N+4);    %N+3..N+4-th columns in PDM matrix is the cost domain
QD=PDM(:,N+5:N+6); %N+5..N+6-th columns in PDM matrix is the quality domain

T=[min(TD,[],2),max(TD,[],2)];   %The first column is the crash/minimal 
                 %durations the last column is the normal/maximal durations
C=[max(CD,[],2),min(CD,[],2)];%The first column is the crash/maximal cost 
                %demands the last column is the normal/minimal cost demands
q=[min(QD,[],2),max(QD,[],2)];%The first column is the crash/minimal 
%quality parameters. The last column is the normal/minimal quality 
%parameters

RD=PDM(:,N+1+3*w:M); 
R=RD; 
for i=0:nR-1     
    R(:,i*w+1:i*w+2)=[max(RD(:,i*w+1:i*w+2),[],2),...
        min(RD(:,i*w+1:i*w+2),[],2)]; 
end                    %Odd columns are the crash/maximal resources demands 
                   %the even columns are the normal/minimal resurce demands
PEM=pem;
P=pem;
Q=ones(N)-P; %Q=1-P

%---End of the initialization of the global values---

n=numel(PEM(PEM>0&PEM<1));  %Number of uncertainties

MinPopSize=100; %Minimal number of population
PopSize=min(max(1000,round((2^n)/100)),MinPopSize); %The population size 
%depends on the number of uncertainties

%---Set of genetic algorithm (GA)---

chromosomec=initpopc(PopSize); %Specify the initial population

% Display:off = There is no need any information to display
% PopulationSize:PopSize = The Population size in a parallel running is
%PopSize
% Vectorized:on and UseParallel:true parallelize the genetic algorithm
% Generations:100 and StallGenLimit:100 sets the maximal number of
%generations to 100
% TolCon:1e-6 means, that the GA will stop, if the given tolerance is
%reached
% CreationFcn:@createpemc specifies createpemc function when generating
%population.
% HybridFcn:@fmincon specifies a hybrid function which is used after
%running GA.
% CrossoverFcn:@crossoverhybridcq specifies a unique crossoverhybridcq 
%crossover function, which recombines two PSM=[DSM,T] matrices  
% MutationFcn:@mutationadaptfeasible randomly generates directions that are
%adaptive with respect to the last successful or unsuccessful generation
% SelectionFcn:{@selectiontournament,4} Tournament selection chooses each 
%parent by choosing Tournament size players at random and then choosing the
%best individual out of that set to be a parent. Tournament size is 4.

options=gaoptimset('Display','off','PopulationSize',PopSize,...
    'Vectorized','on','UseParallel',true,...
    'Generations',100,'StallGenLimit',100,'TolCon',1e-6,...
    'InitialPopulation',chromosomec,'CreationFcn',@createpemc,...
    'HybridFcn',@fmincon,'CrossoverFcn',@crossoverhybridcq,...
    'MutationFcn',@mutationadaptfeasible,...
    'SelectionFcn',{@selectiontournament,4});

LB=[zeros(1,numel(PEM(PEM>0&PEM<1))),min(TD,[],2)']; %Lower Bound
UB=[ones(1,numel(PEM(PEM>0&PEM<1))),max(TD,[],2)'];  %Upper Bound

%---End of setting GA---

%---Runing GA---
% @targetfcncq: the target function (=targetfcncq)
% n+N: the number of variables (n=number of uncertainties + N=number of
%activities
% LB,UB: Lower/Upper Bounds
% @constrainrtparcq: The nonlinear constraint function
% options: The set of options (see above).
%Note: Despite binary and integer values have to be seeked, in order to 
%apply custom creation, crossover and hybrid functions the values
%considered as rounded continuous values

chromosome=ga(@targetfcncq,n+N,[],[],[],[],LB,UB,...
    @constrainrtparcq,options);

%---End of runing GA---  
 
%---Start of output formatting---  

PSM=updatepemc(chromosome); %PSM=[DSM,T] = Logic Domain and the vector 
                            %of the selected modes
                            
PSM(diag(PSM)==0,:)=0; %All excluded tasks dependencies/time demands/ 
PSM(:,diag(PSM)==0)=0; %cost demands are erased from the PSM matrix

PSM=pmtopsmcqr([PSM,zeros(size(PSM,1),1)],T,C,q,R); 
PSM=PSM(:,1:end-1); %cut the last column 

%---Start of Pareto Optimizing---

SST=multires_solve(PSM(:,1:N),PSM(:,N+1),PSM(:,N+4:end)); %Specify an SST 
                    %for a Pareto-Optimal resource allocation 

%---End of Pareto Optimizing---

PSM=[PSM,SST]; %Let SST be the last column

PSM(diag(PSM)==0,:)=0; %All excluded tasks dependencies/time demands/ 
PSM(:,diag(PSM)==0)=0; %cost demands are erased from the PSM matrix

%---End of output formatting---