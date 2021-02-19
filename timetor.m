%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Time to Resource conversion for CTCTP and CTQCTP 
%----------------
%Output:
%R is an N by nR matrix of resource demands
%---------------- 
%Inputs:
%T is an N by 1 vector of task durations
%TD is an N by 2 matrix of min,max time demands
%RD is an N by 2*nR matrix of min,max resource demands
%---------------- 
%Usage:
%R=timetor(T,TD,RD)
%---------------- 
%Example:
%TD=rand(10,2)*20; %Generate Random Time Domain
%RD=rand(10,2*3)*5;%Generate Random Resource Domain
%T=min(TD')'+(max(TD')'-min(TD')').*rand(10,1); %Generate Values Between
%                           %bounds
%R=timetor(T,TD,RD) %Calculate cost demands
%---------------- 
%Prepositions and Requirements:
%1.)TD must be N by 2 RD must be N by 2*nR matrices of non-negative numbers 
%2.)T must be in interval [min(TD')',(max(TD')]
%3.)This model assumes reverse linear bijective function between time and
%resource demands. (If ti<tj=>ri>=rj)

function R=timetor(T,TD,RD)
n=size(RD,2); %n=2*nR;
rmin=[];
rmax=[];
for i=1:2:n-1
    rmin=[rmin,min(RD(:,i:i+1,1)')']; %Select minimal resource demands
    rmax=[rmax,max(RD(:,i:i+1,1)')']; %Select maximal resource demands
end
tmin=repmat(min(TD,[],2),1,1); %Select minimal time demands
tmax=repmat(max(TD,[],2),1,1); %Select maximal time demands
q=tmax-tmin;
R=rmax-rmin;
nR=size(R,2);
for i=1:nR
    R(:,i)=R(:,i).*(tmax-T)./q+rmin(:,i);
end
i=find(q==0);
if numel(i)>0
    R(i)=rmin(i); %If tmax==tmin=>R(i)=rcd min(i)
end