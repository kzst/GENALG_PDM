%Author: Zsolt T. Koszty�n Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Calculate the constraints for a linear continouos time-cost trade-off 
%problem.
%----------------
%Outputs:
%c is a vector of inequities.
%ceq is null vector of Eq.s
%---------------- 
%Input:
%Ce: N by 1 vector of cost demands 
%---------------- 
%Usage:
%[c,ceq]=constrain(Ce)
%---------------- 
%Example:
%This is not executed as a stand-alone fuinction. Instead of specifying
%global values this function is cannot be run.
%---------------- 
%Prepositions and Requirements:
%1.)Global values DSM,T,C,Ct,Cc must be specified.
%2.)DSM should be a binary, upper triangular matrix of the logic domain.

function [c,ceq]=constrain(Ce)
global DSM T C Ct Cc
%global values: 
% DSM is the input logic domain
% T is the column vector of time demands (=the time domain)
% C is the column vector of cost demands (=the cost domain)
% Ct is the time constraint
% Cc is the cost constraint

c(1)=sum(Ce.*floor(diag(DSM)))-Cc; %cost constraints
c(2)=tptfast(DSM,costtotime(Ce,T,C))-Ct; %time constraints
ceq=[];