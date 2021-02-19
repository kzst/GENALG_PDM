%Author: Zsolt T. Koszty�n Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Evaluate EST times of activity
%in order to help the parallelization
%----------------
%Outputs:
%EST: Early Start Time (N by 1 vector)
%----------------
%Inputs:
%DSM: Dependeny Structure Matrix (N by N matrix (logic plan))
%T: Duration time (N by 1 vector)
%----------------
%Usage: 
%EST=est(DSM,T)
%----------------
%Example: - Calculate ESTs
% %    |1,1,0|    |3|
% %DSM=|0,1,0|, T=|4|
% %    |0,0,1|    |5|
%
%EST=est([[1,1,0];[0,1,0];[0,0,1]],[3,4,5]')

function ES=est(DSM,T)
[~,ES]=tptfast(DSM,T);