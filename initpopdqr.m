%Author: Zsolt T. Koszty�n Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Create initial population for resource constrained HDTCTPs and HDTQCTPs
%----------------
%Output:
%chromosome=PopSize by n+2*N matrix of included/excluded dependencies/tasks
%and completion modes and scheduled start times
%---------------- 
%Input:
%PopSize=Size of the population
%---------------- 
%Usage:
%chromosome=initpopdqr(PopSize)
%---------------- 
%Example:
%This is not executed as a stand-alone fuinction. Instead of specifying
%global values this function is cannot be run.
%---------------- 
%Prepositions and Requirements:
%1.)Global values P,T must be specified.

function chromosome=initpopdqr(PopSize)
global P T
%global values: 
% P is the score matrix of inclusion task dependencies/completions 
% T is the column vector of time demands (=the time domain)
nA=sum(diag(P)>0&diag(P)<1); %Number of uncertain acitivieis
n=numel(P(P>0&P<1)); %Number of uncertain activities+uncertain dependencies
N=size(P,1); %Number of activities
w=size(T,2); %Number of modes
chromosomediag=round(rand(PopSize,nA));
chromosome=zeros(PopSize,n+2*N);

for I=1:PopSize

    uncertainties=randoutdiag(P,chromosomediag(I,:)); %First, the set of 
    DSM=updatepem(uncertainties);    %completed activities are specified
    DSM(diag(DSM)==0,:)=0; %Exclude dependencies of excluded tasks
    DSM(:,diag(DSM)==0)=0; %Exclude dependencies of excluded tasks
    MODES=ceil(rand(1,N)*w); %...between in {1,w}
    MODES(diag(DSM)==0)=0; %A duration is set to be zero when a task is 
                           %excluded
    t=zeros(N,1);
    for i=1:N
        if DSM(i,i)==0
            t(i)=0;
        else
            t(i)=T(i,MODES(i));
        end
    end
    [~,EST,~,LST]=tptfast(DSM,t);    
    SST=rand(1,N).*(LST-EST)'+EST'; %Calculation of SST
    chromosome(I,:)=[decision(DSM),MODES,SST]; %chromosome contains
end %uncerain realizations and completion modes and SST