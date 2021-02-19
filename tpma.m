%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Calculate an minimal time demands for a linear continouos time-cost trade-
%off problem.
%----------------
%Output:
%Flag is a scalar value of the optimisation result. In the case of Flag==0,
%t,c are feasible vectors of time, cost demands. If Flag==-1, there is no
%feasible solution.
%t N by 1 vector of minimal time demands.
%c N by 1 vector of cost demands for the minimal time durations.
%---------------- 
%Inputs:
%LD: N by N binary upper triangular matrix of the logic domain 
%---------------- 
%Usage:
%[Flag,t,c]=tpma(LD)
%---------------- 
%Example:
%This is is not executed as a stand-alone fuinction. Instead of specifying
%global values this function is cannot be run.
%---------------- 
%Prepositions and Requirements:
%1.)Global values T,C,Ct,Cc must be specified.
%2.)LD should be a binary, upper triangular matrix of the logic domain.

function [Flag,t,c]=tpma(LD)
global DSM T C Ct Cc
%global values: 
% DSM is the input logic domain
% T is the column vector of time demands (=the time domain)
% C is the column vector of cost demands (=the cost domain)
% Ct is the time constraint
% Cc is the cost constraint

%---Initialization of the global and output varaibles---

Flag=0;
DSM=triu(round(LD));

%---End of the initialization of the global variables---

%---Set of options of fminimax function---
% fiminimiax: fminimax is the applied function to minimize task durations.
% Display:off = There is no need any information to display

options=optimoptions('fminimax','Display','off');

%---End of setting fminimax---

%---Runing fminimax---

c=fminimax(@mpmcost,min(C')',[],[],[],[],min(C')',max(C')',...
    @constrain,options);

% @mpmcost: the target function (=mpmcost)
% min(C')': the x0 vector of initial point
% min(C')',max(C')': linear cost constraints
% @constrain: linear constraint function
% options: The set of options (see above).

%---Evaluation of the results---

t=costtotime(c,T,C);
if sum(c.*floor(diag(DSM)))>Cc+1e-6
    Flag=-1; %There is no optimal solution
end

TPT=tptfast(DSM,t);
if TPT>Ct+1e-6
    Flag=-1; %There is no optimal solution
end