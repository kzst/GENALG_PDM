%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Calculate the constraints for a pareto-optimal resource allocation problem
%----------------
%Outputs:
%c is a vector of inequities.
%ceq is null vector of Eq.s
%---------------- 
%Input:
%SST: N by 1 vector of Scheduled Start Time of tasks 
%---------------- 
%Usage:
%[c,ceq]=constrainrt(SST)
%---------------- 
%Example:
%This is not executed as a stand-alone fuinction. Instead of specifying
%global values this function is cannot be run.
%---------------- 
%Prepositions and Requirements:
%1.)Global values DSM,T,Cr,Ct must be specified.
%2.)DSM should be a binary, upper triangular matrix of the logic domain.

function [c,ceq]=constrainrt(SST)
global DSM T Cr Ct
%global values: 
% DSM is the input logic domain
% T is the column vector of time demands (=the time domain)
% Cr is the resource constraint
% Ct is the time constraint

rMAX=maxres(SST); %Maximal resource constraints
c=zeros(numel(rMAX)+1,1);
TPT=tptsst(DSM,T,SST);
c(1:numel(Cr))=rMAX-Cr; %Maximal resource demands should be lower than Cr.
c(numel(Cr)+1)=TPT-Ct;  %TPT should bo lower than Ct.
ceq=[];