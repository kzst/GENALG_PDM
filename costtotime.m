%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Cost to Time conversion for CTCTP and CTQCTP 
%----------------
%Output:
%T is an N by 1 vector of task durations
%---------------- 
%Inputs:
%C is an N by 1 vector of cost demands
%TD is an N by 2 matrix of min,max time demands
%CD is an N by 2 matrix of min,max cost demands
%---------------- 
%Usage:
%T=costtotime(C,TD,CD)
%---------------- 
%Example:
%TD=rand(10,2)*20; %Generate Random Time Domain
%CD=rand(10,2)*30; %Generate Random Cost Domain
%C=min(CD')'+(max(CD')'-min(CD')').*rand(10,1); %Generate Values Between
%                           %bounds
%T=costtotime(C,TD,CD) %Calculate time demands
%---------------- 
%Prepositions and Requirements:
%1.)TD,CD matrices must be N by 2 matrices of non-negative numbers 
%2.)C must be in interval [min(CD')',(max(CD')]
%3.)This model assumes reverse linear bijective function between time and 
%cost demands. (If ti<tj=>ci>=cj)

function T=costtotime(C,TD,CD)
n=size(C,2); %n=1;
cmin=repmat(min(CD,[],2),1,n); %Select minimal cost demands
cmax=repmat(max(CD,[],2),1,n); %Select maximal cost demands
tmin=repmat(min(TD,[],2),1,n); %Select minimal time demands
tmax=repmat(max(TD,[],2),1,n); %Select maximal time demands
q=cmax-cmin; 
T=(tmax-tmin).*(cmax-C)./q+tmin;
i=find(q==0);
if numel(i)>0
    T(i)=tmin(i); %If cmax==cmin=>T(i)=tmin(i)
end

