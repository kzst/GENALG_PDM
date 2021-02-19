%Author: Zsolt T. Koszty�n Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Implementation of agile project management agent (APMa)
%----------------
%Output:
%PSM: N by M+1 matrix of the calculated project plan. PSM contains the
%logic domain, time and cost demands and (depends on the problem selection 
%it contains) quality and/or resource demands. 
%PSM=[DSM,TD,CD{,QD}{,RD},EST|SST]
% DSM is an N by N upper triangular binary matrix of the calculated logic
% domain
% TD is an N by 1 column vector of task durations 
% CD is an N by 1 column vector of cost demands
%depends on the problem selection
% QD is an N by 1 column vector of quality demands
% RD is an N by nR matrix of resource demands
%depends on the resource allocation is performed or not
% SST is an N by 1 column vector of scheduled start time of tasks if
%resource allocation is applied
% EST is an N by 1 column vector of earliest start time of tasks if
%resource demands are not considered
%---------------- 
%Inputs:
%PDM: N by N+(M-N)*W matrix of the stochastic project plan. 
%PDM=[PEM,TD,CD{,QD}{,RD}], where PEM is an N by N upper triangular matrix 
%of logic domain, 
% TD is an N by w matrix of task durations,
% CD is an N by w matrix of cost demands,
%depends on the problem selection
% QD is an N by w matrix of quality parameters,
%depends on the resource allocation is performed or not
% RD is an N by w*nR matrix of resource demands.
%const: 3+nR+1 elements row vector of [Ct,Cc,Cq,CR,Cs] values, where
% Ct is the time constraint (max constraint)
% Cc is the cost constraint (max constraint)
% Cq is the quality constraint (min constraint)
% CR are the resource constraints (max constraint)
% Cs is the score constraint (min constraint)
%CONST: 3..3+nR+1 (depends on the problem selection) elements row vector of 
%[Ct,Cc,{Cq,CR,}Cs] values, where
% Ct is the time constraint (max constraint, required in all problems)
% Cc is the cost constraint (max constraint, required in all problems)
% Cq is the quality constraint (min constraint, required in TQCTP problems))
% CR are the resource constraints (max constraint, required in 
%resource-allocation problems)
% Cs is the score constraint (min constraint, required in all problems)
%Select: 1-18 scalar of problem selection
% 1,3,4: APSP: Agile Project Scheduling Problem
% 2,5,6: APSPq: Agile Project Scheduling Problem with quality parameters
% 7,9,10: RC-APSP: Resource-Constraint Agile Project Scheduling Problem 
% 8,11,12: PO-APSP: Pareto-Optimal resource-allocated Agile Project 
%Scheduling Problem 
%13,15,17: RC-APSPq: Resource-Constraint Agile Project Scheduling Problem 
%with quality parameters
%14,16,18: PO-APSPq: Pareto-Optimal resource-allocated Agile Project Scheduling 
%with quality parameters
%typetfcn: Type of target function 
% {0=maxTPQ,} 1=minTPT, 2=minTPC, 3=maxTPS, {4=minUF,} ~ composite
%w: Number of modes
%---------------- 
%Usage:
%PSM=mapma(PDM,CONS,Select,typefcn,w)
%---------------- 
%Example:
%PDM=[triu(rand(10)*.5+.5),20*rand(10,3),30*rand(10,3),rand(10,3),...
%   5*rand(10,6)];
%const=[percentt(PDM,3,.9),percentc(PDM,3,.9),percentqr(PDM,3,1)',...
%   percentq(PDM,3,.7),percents(PDM,.7)];
%Select=15;
%typefcn=999; %let target function be a composite target function
%w=3;
%tic;PSM=mapma(PDM,CONS,Select,typefcn,w);toc
%---------------- 
%Prepositions and Requirements:
%1.)Agile Project Management agent does not consider completion modes. If
%there are more than one mode, a normal duration/cost/quality is
%considered.
%2.)The logic domian of the input PDM matrix must be an upper triangular 
%matrix, where the matrix elements are between 0 to 1 interval.
%3.)At least one matrix element is lower than 1 and upper than 0. => There
%are at least one decision variable
%4.)The number of resources (nR) is a positive number
%5.)The elements of Time/Cost/Quality/Resource Domains are positive real 
%numbers.

function PSM=mapma(PDM,CONS,Select,typefcn,w)
N=size(PDM,1); %Number of activities
PEM=PDM(:,1:N); %PEM is the first N by N domain of logic plan
M=size(PDM,2); %M is the number of activities + modes of activities, costs, 
%(quality parameters, resources (depends on the problem))
switch Select
    case {1,3,4} % No resources, no alternatives/select normal parameters, 
        %no quality parameters; select variables
        %Despite w=1 is assumed, this section works, if w>1, however, in
        %this case normal paramteres (i.e. duration, cost demands) are
        %selected to evaluate
        [TD,MODES]=max(PDM(:,N+1:N+w)'); %Normal=maximal duration is 
        %selected.
        TD=TD';
        C=PDM(:,N+w+1:N+2*w);
        CD=zeros(N,1);
        MODES=MODES';
        for I=1:N %The cost of normal duration is the normal cost
            if MODES(I)>0
                CD(I)=C(I,MODES(I));
            end
        end
        PDM=[PEM,TD,CD];
        PSM=apmgen(PDM,CONS,typefcn);
   case {2,5,6} % No resources, no alternatives/select normal parameters, 
        %with quality parameters; select variables
        %Despite w=1 is assumed, this section works, if w>1, however, in
        %this case normal paramteres (i.e. duration, cost demands, quality 
        %parameters) are selected to evaluate
        [TD,MODES]=max(PDM(:,N+1:N+w)'); %Normal=maximal duration is 
        %selected.
        TD=TD';
        C=PDM(:,N+w+1:N+2*w);
        Q=PDM(:,N+2*w+1:N+3*w);
        CD=zeros(N,1);
        QD=zeros(N,1);
        MODES=MODES';
        for I=1:N
            if MODES(I)>0
                CD(I)=C(I,MODES(I));
                QD(I)=Q(I,MODES(I));
            end
        end % The cost of normal duration is the normal cost. Similarly the
        %the normal quality parameteds is the quality paremeters of normal
        %durations
        PDM=[PEM,TD,CD,QD];
        PSM=apmgenq(PDM,CONS,typefcn);
    case {7,9,10} % With resources, select normal parameters, no quality 
        %parameters; select variables
        %Despite w=1 is assumed, this section works, if w>1, however, in
        %this case normal paramteres (i.e. duration, cost demands, quality 
        %parameters) are selected to evaluate
        [TD,MODES]=max(PDM(:,N+1:N+w)'); %Normal=maximal duration is 
        %selected.
        TD=TD';
        C=PDM(:,N+w+1:N+2*w);
        R=PDM(:,N+2*w+1:end);
        nR=(M-N-2*w)/w;
        CD=zeros(N,1);
        RD=zeros(N,nR);
        MODES=MODES';
        for I=1:N % The cost/resource demands of normal duration is the 
            %normal cost/resources
            if MODES(I)>0
                CD(I)=C(I,MODES(I));
                for J=1:nR
                    RD(I,J)=R(I,(J-1)*MODES(I)+1);
                end
            end
        end
        PDM=[PEM,TD,CD,RD];
        PSM=hpmgendr(PDM,CONS,typefcn);     
    case {8,11,12} % Pareto-optimal resource allocation, select normal 
        %parameters, no quality parameters; select variables
        %Despite w=1 is assumed, this section works, if w>1, however, in
        %this case normal paramteres (i.e. duration, cost demands, quality 
        %parameters) are selected to evaluate
        [TD,MODES]=max(PDM(:,N+1:N+w)'); %Normal=maximal duration is 
        %selected.
        TD=TD';
        C=PDM(:,N+w+1:N+2*w);
        R=PDM(:,N+2*w+1:end);
        nR=(M-N-2*w)/w;
        CD=zeros(N,1);
        RD=zeros(N,nR);
        MODES=MODES';
        for I=1:N %The cost/resource demands of normal duration is the
            %normal cost/resources
            if MODES(I)>0
                CD(I)=C(I,MODES(I));
                for J=1:nR
                    RD(I,J)=R(I,(J-1)*MODES(I)+1);
                end
            end
        end
        PDM=[PEM,TD,CD,RD];
        PSM=apmgenpr(PDM,CONS,typefcn);  
    case {13,15,17} % With resources, select normal parameters, with quality parameters;
        %select variables
        %Despite w=1 is assumed, this section works, if w>1, however, in
        %this case normal paramteres (i.e. duration, cost demands, quality 
        %parameters, resource demands) are selected to evaluate
        [TD,MODES]=max(PDM(:,N+1:N+w)'); %Normal=maximal duration is 
        %selected.
        TD=TD';
        C=PDM(:,N+w+1:N+2*w);
        Q=PDM(:,N+2*w+1:N+3*w);
        R=PDM(:,N+3*w+1:end);
        nR=(M-N-3*w)/w;
        CD=zeros(N,1);
        QD=zeros(N,1);
        RD=zeros(N,nR);
        MODES=MODES';
        for I=1:N %The cost/resource demands of normal duration is the
            %normal cost/resources
            if MODES(I)>0
                CD(I)=C(I,MODES(I));
                QD(I)=Q(I,MODES(I));
                for J=1:nR
                    RD(I,J)=R(I,(J-1)*MODES(I)+1);
                end
            end
        end
        PDM=[PEM,TD,CD,QD,RD];
        PSM=hpmgendqr(PDM,CONS,typefcn); 
    case {14,16,18} % Pareto-optimal resource allocation, select normal 
        %parameters, with quality parameters; select variables
        %Despite w=1 is assumed, this section works, if w>1, however, in
        %this case normal paramteres (i.e. duration, cost demands, quality 
        %parameters, resource demands) are selected to evaluate
        [TD,MODES]=max(PDM(:,N+1:N+w)'); %Normal=maximal duration is 
        %selected.
        TD=TD';
        C=PDM(:,N+w+1:N+2*w);
        Q=PDM(:,N+2*w+1:N+3*w);
        R=PDM(:,N+3*w+1:end);
        nR=(M-N-3*w)/w;
        CD=zeros(N,1);
        QD=zeros(N,1);
        RD=zeros(N,nR);
        MODES=MODES';
        for I=1:N %The cost/resource demands of normal duration is the
            %normal cost/resources
            if MODES(I)>0
                CD(I)=C(I,MODES(I));
                QD(I)=Q(I,MODES(I));
                for J=1:nR
                    RD(I,J)=R(I,(J-1)*MODES(I)+1);
                end
            end
        end
        PDM=[PEM,TD,CD,QD,RD];
        PSM=apmgenqpr(PDM,CONS,typefcn); 
    case {19,20} % With resources, select normal parameters, no quality 
        %parameters; select variables
        %Despite w=1 is assumed, this section works, if w>1, however, in
        %this case normal paramteres (i.e. duration, cost demands, quality 
        %parameters) are selected to evaluate
        [TD,MODES]=max(PDM(:,N+1:N+w)'); %Normal=maximal duration is 
        %selected.
        TD=TD';
        C=PDM(:,N+w+1:N+2*w);
        R=PDM(:,N+2*w+1:end);
        nR=(M-N-2*w)/w;
        CD=zeros(N,1);
        RD=zeros(N,nR);
        MODES=MODES';
        for I=1:N % The cost/resource demands of normal duration is the 
            %normal cost/resources
            if MODES(I)>0
                CD(I)=C(I,MODES(I));
                for J=1:nR
                    RD(I,J)=R(I,(J-1)*MODES(I)+1);
                end
            end
        end
        PDM=[PEM,TD,CD,RD];
        PSM=parhpmgendr(PDM,CONS);
    case {21,22} %With resources, select normal parameters, with quality parameters;
        %select variables
        %Despite w=1 is assumed, this section works, if w>1, however, in
        %this case normal paramteres (i.e. duration, cost demands, quality 
        %parameters, resource demands) are selected to evaluate
        [TD,MODES]=max(PDM(:,N+1:N+w)'); %Normal=maximal duration is 
        %selected.
        TD=TD';
        C=PDM(:,N+w+1:N+2*w);
        Q=PDM(:,N+2*w+1:N+3*w);
        R=PDM(:,N+3*w+1:end);
        nR=(M-N-3*w)/w;
        CD=zeros(N,1);
        QD=zeros(N,1);
        RD=zeros(N,nR);
        MODES=MODES';
        for I=1:N %The cost/resource demands of normal duration is the
            %normal cost/resources
            if MODES(I)>0
                CD(I)=C(I,MODES(I));
                QD(I)=Q(I,MODES(I));
                for J=1:nR
                    RD(I,J)=R(I,(J-1)*MODES(I)+1);
                end
            end
        end
        PDM=[PEM,TD,CD,QD,RD];
        PSM=parhpmgendqr(PDM,CONS);    
    otherwise
        PSM=PDM;
end
