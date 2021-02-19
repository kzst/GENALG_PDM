%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Simulating the effects of the changes of customer claims
%Two steps of Monte-Carlo simulation is applied for selecting (1) and
%changing (2) task demands and the structures.
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
%PDMout=phase3(PDM,p,s)
%---------------- 
%Example:
%PDM=generatepdmq(30,0.05,0,20,30,20,2,2,2);
%a=-.1;%New values will be greater or equal than 90% of the original values
%b=0.3;%New values will be lower or equal than 90% of the original values
%PDM1=phase1(PDM,a,b)
%p=.05;%Probability factor for task selection
%s=2.0;%Scale factor: the ratio of the modification
%PDM2=phase2(PDM1,p,s)
%PDM3=phase2(PDM2,p,s)
%---------------- 
%Prepositions and Requirements:
%1.)All domains can be changed

function PDMout=phase3_restricted(PDM,p,s,N)
n=size(PDM,1);
m=size(PDM,2);
PDMout=PDM;
PDMout(1:n,1:n)=max(min(ones(n),PDM(1:n,1:n)+...
    s*triu((rand(n)<p).*rand(n))),zeros(n)); %New dependencies/task 
if m>n  %occurances is generated with probability value p
    Z=zeros(n,m-n);
    PDMout(:,n+1:m)=PDMout(:,n+1:m)+(PDM(:,n+1:m)==Z).*rand(n,m-n,...
        'like',PDM(diag(PDM)~=0,n+1:m)).*... %If new task occurances 
        repmat(mean(PDM(diag(PDM)~=0,n+1:m)),n,1); %is generated then
    for i=1:n                %demands will be similar to the other demands
        for j=n+1:m
            if ((PDM(i,j)>0)&&(PDM(i,j)<=1)&&(PDMout(i,j)>1))
                PDMout(i,j)=1; %Quality should not be greater than 1
            end
        end
    end
end
PDMout(1:N,1:N)=PDM(1:N,1:N);
PDMout(diag(PDMout)==0,:)=0;   %Exluded task demands are also excluded
PDMout(1:n,diag(PDMout)==0)=0; %Exluded task demands are also excluded
