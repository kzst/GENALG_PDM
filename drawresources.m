%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Draw the resource diagram
%----------------
%Input:
%DSM: N by N upper triangular binary matrix of logic domain
%T: N by 1 vector of time demands
%R: N by nR matrix of resource demands
%SST: N by 1 vector of scheduled start time
%---------------- 
%Usage:
%drawresources(DSM,T,R,SST)
%---------------- 
%Example: draw resource demands for tasks
%DSM=triu(round(rand(10))); T=rand(10,1)*20; R=rand(10,3)*5; %generate
%initial domains
%[~,EST,~,LST,~]=tptfast(DSM,T); %calculate EST and LST
%SST=EST+rand(10,1).*(LST-EST); %generate SST between EST and LST
%drawresources(DSM,T,R,SST); %draw resource graphs

function drawresources(DSM,T,R,SST)
[~,EST,~,LST,~]=tptfast(DSM,T);
[BP,RESFUNC]=resfunc(DSM,EST,T,R);
figure('Name','EST');
drawres(BP,RESFUNC);
[BP,RESFUNC]=resfunc(DSM,LST,T,R);
figure('Name','LST');
drawres(BP,RESFUNC);
[BP,RESFUNC]=resfunc(DSM,SST,T,R);
figure('Name','SST');
drawres(BP,RESFUNC);
