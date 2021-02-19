%Author: Zsolt T. Kosztyan Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Simulating the effects of the uncertainty of specification
%Random modification, which follows beta distribution for 
%time/cost/resource demands and quality parameters
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
%sched=ratio of scheduled tasks
%---------------- 
%Usage:
%PDMout=phase1betasched(PDM,a,b,sched)
%---------------- 
%Example:
%PDM=generatepdmq(30,0.05,0,20,30,20,2,2,2);
%a=-.1;%New values will be greater or equal than 90% of the original values
%b=0.3;%New values will be lower or equal than 90% of the original values
%PDM1b=phase1betasched(PDM,a,b,sched)
%---------------- 
%Prepositions and Requirements:
%1.)Only TD,CD,QD and RD domains will be modified.

function PDMout=phase1betasched(PDM,a,b,sched)
n=size(PDM,1);
t=(a+b)/6;
r1=a/t;
r2=b/t;
pPDM=[PDM(1+round((n-1)*sched):n,1+round((n-1)*sched):end)];
if size(pPDM,1)>1
    PDMout=pPDM;
    np=size(pPDM,1);
    m=size(pPDM,2);
    if b>a
        alpha=6*(1-r1)/(r2-r1);
        beta=6*(r2-1)/(r2-r1);
        M=b-a;
        if m>np
            PDMout(:,np+1:end)=pPDM(:,np+1:end)+(M*random('beta',alpha,beta,...
                np,m-np)+a).*pPDM(:,np+1:end);
            for i=1:np
                for j=np+1:m
                    if ((pPDM(i,j)>0)&&(pPDM(i,j)<=1)&&(PDMout(i,j)>1))
                        PDMout(i,j)=1; %Quality should not be greater than 1
                    end
                end
            end
        end
    end
    temp=PDMout;
    PDMout=PDM;
    PDMout(1+round((n-1)*sched):n,1+round((n-1)*sched):end)=temp;
else
    PDMout=PDM;
end