%Author: Zsolt T. Koszty�n Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Calculate an optimal project structure matrix and optimal resource
%allocation for a resource constraint hybrid continouos time-cost trade-off
%problem (RCHCTCTP)
%The algorithm is the implementation of resource-constrained Hybrid Project 
%Management agent (rcHPMa) by hybrid genetic algorithm 
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
%PDM: PDM=[PEM,TD,CD,RD] matrix (N by N+(2+nR)*2), 
%where PEM is an N by N upper triangular matrix of logic domain, 
% TD is an N by 2 matrix of task durations
% CD is an N by 2 matrix of cost demands
% RD is an N by 2*nR matrix of resource demands
%const: 2+nR+1 elements row vector of [Ct,Cc,CR,Cs] values, where
% Ct is the time constraint (max constraint)
% Cc is the cost constraint (max constraint)
% CR are the resource constraints (max constraint)
% Cs is the score constraint (min constraint)
%typefcn: Type of target function 
% 1=minTPT, 2=minTPC, 3=maxTPS, 4=minUF, ~ composite
%---------------- 
%Usage:
%PSM=hpmgencr(PDM,const,typefcn)
%---------------- 
%Example:
%PDM=[triu(rand(10)*.5+.5),20*rand(10,2),30*rand(10,2),5*rand(10,6)];
%const=[percentt(PDM,2,.9),percentc(PDM,2,.9),percentr(PDM,2,1)',...
%   percents(PDM,.7)];
%typefcn=999; %let target function be a composite target function
%tic;PSM=hpmgencr(PDM,const,typefcn);toc
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
%modes. Therefore the algorithm match maximal cost/resource demands to the 
%minimal durations and minimal cost/resource demands to themaximal 
%durations. 

function PSM=hpmgencr(PDM,const,typefcn)
global PEM P Q T C R Ct Cc CR Cs Typ

%global values: 
% PEM is the input logic domain
% P is the score matrix of inclusion task dependencies/completions 
%the logic domain (initially PEM=P)
% Q is the score matrix of exclusion task dependencies/completions. 
%In this simulation Q=1-P
% T is the N by 2 matrix of time demands (=the time domain)
% C is the N by 2 matrix of cost demands (=the cost domain)
% R is the N by 2*nR matrix of resource demands (=the resource domain)
% Ct is the time constraint
% Cc is the cost constraint
% CR is the row vector of resource constraints
% Cs is the score constraint
% Typ is the selection number of the target functions (see above)

%---Initialization of the global values---

w=2; %In the case of CTCTP, the number of modes(w) is always 2.

Typ=typefcn;
Ct=const(1);
Cc=const(2);
CR=const(3:end-1);
Cs=const(end);
N=size(PDM,1); %The number of activities 
M=size(PDM,2); %M=N+2*2+2*nR
nR=(M-N-4)/2; %nR: the number of resources 
pem=triu(PDM(:,1:N)); %Only the upper triangle will be considered
TD=PDM(:,N+1:N+2);    %N+1..N+2-th columns in PDM matrix is the time domain
CD=PDM(:,N+3:N+4);    %N+3..N+4-th columns in PDM matrix is the cost domain
RD=PDM(:,N+1+2*w:end); %N+5..M-th columns in PDM matrix is the resource 
                       %domain
T=[reshape(min(TD,[],2),[],1),reshape(max(TD,[],2),[],1)];   %The first column is the crash/minimal 
               %durations the last column is the normal/maximal durations
C=[reshape(max(CD,[],2),[],1),reshape(min(CD,[],2),[],1)];%The first column is the crash/maximal cost
               % demands the last column is the normal/minimal cost demands

R=RD;
for i=0:nR-1
    R(:,i*w+1:i*w+2)=[max(RD(:,i*w+1:i*w+2)')',min(RD(:,i*w+1:i*w+2)')'];
end                    %Odd columns are the crash/maximal resources demands 
                   %the even columns are the normal/minimal resurce demands
PEM=pem;
P=pem;
Q=ones(N)-P; %Q=1-P

%---End of the initialization of the global values---

n=numel(PEM(PEM>0&PEM<1));  %Number of uncertainties

MinPopSize=100; %Minimal number of population
PopSize=min(max(1000,round((2^n)/100)),MinPopSize);    %The population size 
%depends on the number of uncertainties

%---Set of genetic algorithm (GA)---
chromosomecr=initpopcr(PopSize); %Specify the initial population

% Display:off = There is no need any information to display
% PopulationSize:PopSize = The Population size in a parallel running is
%PopSize
% Vectorized:on and UseParallel:true parallelize the genetic algorithm
% Generations:100 and StallGenLimit:100 sets the maximal number of
%generations to 100
% TolCon:1e-6 means, that the GA will stop, if the given tolerance is
%reached
% CreationFcn:@createpemcr specifies createpemcr function when generating
%population.
% HybridFcn:@fmincon specifies a hybrid function which is used after
%running GA.
% CrossoverFcn:@crossoverhybridcr specifies a unique crossoverhybridcr 
%crossover function, which recombines two PSM=[DSM,T,SST] matrices  
% MutationFcn:@mutationadaptfeasible randomly generates directions that are
%adaptive with respect to the last successful or unsuccessful generation
% SelectionFcn:{@selectiontournament,4} Tournament selection chooses each 
%parent by choosing Tournament size players at random and then choosing the
%best individual out of that set to be a parent.�Tournament size�is�4.

options=gaoptimset('Display','off','PopulationSize',PopSize,...
    'Vectorized','on','UseParallel',true,...
    'Generations',100,'StallGenLimit',100,'TolCon',1e-6,...
    'InitialPopulation',chromosomecr,'CreationFcn',@createpemcr,...
    'HybridFcn',@fmincon,'CrossoverFcn',@crossoverhybridcr,...
    'MutationFcn',@mutationadaptfeasible,...
    'SelectionFcn',{@selectiontournament,4});

%---Setting the global lower/upper bounds---

[~,EST]=tptfast(floor(PEM),min(T,[],2)'); %EST when all uncertainties are 
                                      %excluded from the DSM
[~,~,~,LST]=tptfast(ceil(PEM),max(T,[],2)');%LST when all uncertainties are 
                                         %completed

LST(LST<0)=0;
Z=numel(PEM(PEM>0&PEM<1));

if Z>0
    LB=[zeros(1,Z),reshape(min(TD,[],2),1,[]),EST'];      %Lower Bound
    UB=[ones(1,Z),reshape(max(TD,[],2),1,[]),LST'];        %Upper Bound
else
    LB=[reshape(min(TD,[],2),1,[]),EST'];      %Lower Bound
    UB=[reshape(max(TD,[],2),1,[]),LST'];      %Upper Bound

end
%---End of setting GA---

%---Runing GA---
% @targetfcncr: the target function (=targetfcncr)
% n+2*N: the number of variables (n=number of uncertainties + N=number of
%activities
% LB,UB: Lower/Upper Bounds
% @constrainrtparcr: The nonlinear constraint function
% options: The set of options (see above).
%Note: Despite binary and integer values have to be seeked, in order to 
%apply custom creation, crossover and hybrid functions the values
%considered as rounded continuous values

chromosome=ga(@targetfcncr,n+2*N,[],[],[],[],LB,UB,...
    @constrainrtparcr,options);

%---End of runing GA---  
 
%---Start of output formatting---  

PSM=updatepemcr(chromosome); %PSM=[DSM,T,SST] = Logic Domain and the vector 
                            %of the selected modes and the vector of
                            %scheduled start time of activities
                            
PSM(diag(PSM)==0,:)=0; %All excluded tasks dependencies/time demands/ 
PSM(:,diag(PSM)==0)=0; %cost demands are erased from the PSM matrix

PSM=pmtopsmcr(PSM,T,C,R); %The output PSM = [DSM,TD,CD,RD,SST]

PSM(diag(PSM)==0,:)=0; %All excluded tasks dependencies/time demands/ 
PSM(:,diag(PSM)==0)=0; %cost demands are erased from the PSM matrix

%---End of output formatting---