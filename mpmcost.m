%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Target function of TCTP: Minimize durations
%----------------
%Output:
%t: is an N by 1 vector of task duration
%---------------- 
%Inputs:
%c: N by 1 vector of cost demands
%---------------- 
%Usage:
%t=mpmcost(c)
%---------------- 
%Example:
%This is is not executed as a stand-alone fuinction. Instead of specifying
%global values this function is cannot be run.
%---------------- 
%Prepositions and Requirements:
%1.)Global values T,C,Ct,Cc must be specified.
%2.)LD should be a binary, upper triangular matrix of the logic domain.

function t=mpmcost(c)
global DSM T C
%global values: 
% DSM is the input logic domain
% T is the column vector of time demands (=the time domain)
% C is the column vector of cost demands (=the cost domain)

t=tptfast(DSM,costtotime(c,T,C));