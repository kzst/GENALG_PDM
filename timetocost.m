%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Time to Cost conversion for CTCTP and CTQCTP 
%----------------
%Output:
%C is an n by 1 vector of cost demands
%---------------- 
%Inputs:
%T is an n by 1 vector of task durations
%TD is an n by 2 matrix of min,max time demands
%CD is an n by 2 matrix of min,max cost demands
%---------------- 
%Usage:
%C=timetocost(T,TD,CD)
%---------------- 
%Example:
%TD=rand(10,2)*20; %Generate Random Time Domain
%CD=rand(10,2)*30; %Generate Random Cost Domain
%T=min(TD')'+(max(TD')'-min(TD')').*rand(10,1); %Generate Values Between
%                           %bounds
%C=timetocost(T,TD,CD) %Calculate cost demands
%---------------- 
%Prepositions and Requirements:
%1.)TD,CD matrices must be N by 2 matrices of non-negative numbers 
%2.)T must be in interval [min(TD')',(max(TD')]
%3.)This model assumes reverse linear bijective function between time and
%cost demands. (If ti<tj=>ci>=cj)

function C=timetocost(T,TD,CD)
n=size(T,2); %n=1
cmin=repmat(min(CD,[],2),1,n); %Select minimal cost demands
cmax=repmat(max(CD,[],2),1,n); %Select maximal cost demands
tmin=repmat(min(TD,[],2),1,n); %Select minimal time demands
tmax=repmat(max(TD,[],2),1,n); %Select maximal time demands
q=tmax-tmin;
C=(cmax-cmin).*(tmax-T)./q+cmin;
i=find(q==0);
if numel(i)>0
    C(i)=cmin(i); %If tmax==tmin=>C(i)=cmin(i)
end

