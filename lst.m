%Author: Zsolt T. Kosztyan Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Evaluate LST times of activity
%in order to help the parallelization
%----------------
%Outputs:
%LST: Late Start Time (N by 1 vector)
%----------------
%Inputs:
%DSM: Dependeny Structure Matrix (N by N matrix (logic plan))
%T: Duration time (N by 1 vector)
%----------------
%Usage: 
%LST=lst(DSM,T)
%----------------
%Example: - Calculate LSTs
% %    |1,1,0|    |3|
% %DSM=|0,1,0|, T=|4|
% %    |0,0,1|    |5|
%
%LST=lst([[1,1,0];[0,1,0];[0,0,1]],[3,4,5]')

function LS=lst(DSM,T)
[~,~,~,LS]=tptfast(DSM,T);