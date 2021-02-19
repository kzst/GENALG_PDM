%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Time to Quality conversion for CTQCTP 
%----------------
%Output:
%Q is an n by 1 vector of quality paramters
%---------------- 
%Inputs:
%T is an n by 1 vector of task durations
%TD is an n by 2 matrix of min,max time demands
%QD is an n by 2 matrix of min,max quality parameters
%---------------- 
%Usage:
%Q=timetoquality(T,TD,QD)
%---------------- 
%Example:
%TD=rand(10,2)*20; %Generate Random Time Domain
%QD=rand(10,2); %Generate Random Cost Domain
%T=min(TD')'+(max(TD')'-min(TD')').*rand(10,1); %Generate Values Between
%                           %bounds
%Q=timetoquality(T,TD,QD) %Calculate cost demands
%---------------- 
%Prepositions and Requirements:
%1.)TD,QD matrices must be N by 2 matrices of non-negative numbers 
%2.)T must be in interval [min(TD')',(max(TD')]
%3.)This model assumes linear bijective function between time and
%quality demands. (If ti<tj=>qi<=qj)

function Q=timetoquality(T,TD,QD)
N=size(T,2);
qmin=repmat(min(QD,[],2),1,N); %Select minimal quality parameters
qmax=repmat(max(QD,[],2),1,N); %Select maximal quality parameters
tmin=repmat(min(TD,[],2),1,N); %Select minimal time demands
tmax=repmat(max(TD,[],2),1,N); %Select maximal time demands
q=tmax-tmin;
Q=(qmax-qmin).*(T-tmin)./q+qmin; 
i=find(q==0);
if numel(i)>0
    Q(i)=qmin(i); %If tmax==tmin=>Q(i)=qmin(i)
end

