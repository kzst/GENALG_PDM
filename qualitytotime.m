%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Quality to Time conversion for CTQCTP 
%----------------
%Output:
%T is an N by 1 vector of task durations
%---------------- 
%Inputs:
%Q is an N by 1 vector of quality parameters
%TD is an N by 2 matrix of min,max time demands
%QD is an N by 2 matrix of min,max quality parameters
%---------------- 
%Usage:
%T=qualitytotime(Q,TD,QD)
%---------------- 
%Example:
%TD=rand(10,2)*20; %Generate Random Time Domain
%QD=rand(10,2); %Generate Random Cost Domain
%Q=min(QD')'+(max(QD')'-min(QD')').*rand(10,1); %Generate Values Between
%                           %bounds
%T=qualitytotime(Q,TD,QD) %Calculate cost demands
%---------------- 
%Prepositions and Requirements:
%1.)TD,QD matrices must be N by 2 matrices of non-negative numbers 
%2.)T must be in interval [min(TD')',(max(TD')]
%3.)This model assumes linear bijective function between time and
%quality demands. (If ti<tj=>qi<=qj)

function T=qualitytotime(Q,TD,QD)
N=size(Q,2);
qmin=repmat(min(QD,[],2),1,N); %Select minimal quality parameters
qmax=repmat(max(QD,[],2),1,N); %Select maximal quality parameters
tmin=repmat(min(TD,[],2),1,N); %Select minimal time demands
tmax=repmat(max(TD,[],2),1,N); %Select maximal time demands
q=qmax-qmin;
T=(tmax-tmin).*(Q-qmin)./q+tmin;
i=find(q==0);
if numel(i)>0
    T(i)=tmin(i); %If qmax==qmin=>T(i)=tmin(i)
end
