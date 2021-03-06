%Author: Zsolt T. Koszty�n Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Simulating the effects of the extreme changes of demands
%Two steps of Monte-Carlo simulation is applied for selecting (1) and
%changing (2) task demands
%----------------
%Output:
%PDMout: PDM matrix with same structure as the input PDM matrix
%---------------- 
%Input:
%PDM=[PEM,TD,CD,{QD,RD}]: Project Domain Matrix
%where PEM is an N by N upper triangular matrix of logic domain, 
% TD is an N by w matrix of task durations
% CD is an N by w matrix of cost demands
% {QD is an N by w matrix of quality parameters}} %optional
% {RD is an N by w*nR matrix of resource demands} %optional
%a=ratio of minimal value
%b=ratio of maximal value
%---------------- 
%Usage:
%PDMout=phase2(PDM,p,s)
%---------------- 
%Example:
%PDM=generatepdmq(30,0.05,0,20,30,20,2,2,2);
%a=-.1;%New values will be greater or equal than 90% of the original values
%b=0.3;%New values will be lower or equal than 90% of the original values
%PDM1=phase1(PDM,a,b)
%p=.05;%Probability factor for task selection
%s=2.0;%Scale factor: the ratio of the modification
%PDM2=phase2(PDM1,p,s)
%---------------- 
%Prepositions and Requirements:
%1.)Only TD,CD,QD and RD domains will be modified.

function PDMout=phase2(PDM,p,s)
n=size(PDM,1);
m=size(PDM,2);
PDMout=PDM;
if m>n
    PDMout(:,n+1:end)=PDM(:,n+1:end)+(rand(n,m-n)<p).*(PDM(:,n+1:end)*s);
    for i=1:n
        for j=n+1:m
            if ((PDM(i,j)>0)&&(PDM(i,j)<=1)&&(PDMout(i,j)>1))
                PDMout(i,j)=1; %Quality should not be greater than 1
            end
        end
    end
end