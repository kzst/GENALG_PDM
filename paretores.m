%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Calculate the Pareto-Optimal Resource Allocation.
%----------------
%Output:
%SST N by 1 vector of Scheduled Start Time of tasks.
%---------------- 
%Inputs:
%LD: N by N binary upper triangular matrix of the logic domain 
%TD: N by 1 vector of the time domain 
%RD: N by nR matrix of the resource demands 
%---------------- 
%Usage:
%SST=paretores(LD,TD,RD)
%---------------- 
%Example:
%LD=triu(round(rand(10))); %Specify Logic Domain
%TD=rand(10,1)*30; %Specify Time Domain
%RD=rand(10,3)*5; %Specify Resource Domain
%tic;SST=paretores(LD,TD,RD);toc
%---------------- 
%Prepositions and Requirements:
%1.)LD should be a binary, upper triangular matrix of the logic domain.
%2.)TD,RD contains positive values.

function SST=paretores(LD,TD,RD)
global DSM T R
%global values: 
% DSM is the input logic domain
% T is the column vector of time demands (=the time domain)
% R is the matrix of resource demands (=the resource domain)

%---Initialization of the global varaibles---

DSM=triu(round(LD)); %DSM must be an upper triangular binary matrix
T=TD;
R=RD;
N=numel(TD);

%---End of the initialization of the global variables---

[~,EST,~,LST,~]=tptfast(LD,TD); %Calculating EST and LST, 
                                %where EST<=SST<=LST
                                
%---Set of multi-objective genetic algorithm (MOGA)---

% Vectorized:off and UseParallel:false serialize the genetic algorithm
% Display:off = There is no need any information to display
                                
options=gaoptimset('Vectorized','off','UseParallel',false,'Display','off');

%---End of setting MOGA---

%---Runing MOGA---
% @maxres: the target function (=maxres, minimize the maximal resource 
%demands)
% N: the number of variables
% EST,LST: lower, upper bounds. 
% options: The set of options (see above).

X=gamultiobj(@maxres,N,[],[],[],[],EST,LST,options);

%---End of runing MOGA---

%---Start of output formatting---

if numel(X)>0
    SST=X(1,:)'; %SST will be the first pareto-optimal solution.
else
    SST=EST; %SST=EST if there is no optimal solution.
end