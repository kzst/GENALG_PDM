%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Draw the resource diagram
%----------------
%Input:
%BP: vector of breakpoints in resource diagram
%RESFUNC: matrix of resource demands (resource demands in breakpoints)
%---------------- 
%Usage:
%drawres(BP,RESFUNC)
%---------------- 
%Example: draw resource demands for tasks, which scheduled as early as
%possible
%DSM=triu(round(rand(10))); T=rand(10,1)*20; R=rand(10,3)*5; %generate
%initial domains
%[~,EST,~,~,~]=tptfast(DSM,T); %calculate EST
%[BP,RESFUNC]=resfunc(DSM,EST,T,R); %calculate BP and RESFUNC
%figure('Name','EST');
%drawres(BP,RESFUNC);

function drawres(BP,RESFUNC)
resources=numel(RESFUNC(1,:));
for i=1:resources
    subplot(resources,1,i)
    stairs(BP,RESFUNC(:,i));
    title(strcat('r_',num2str(i)));
    xlabel('day');
    ylabel('resource');
end
end