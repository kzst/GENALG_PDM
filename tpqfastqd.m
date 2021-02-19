%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Calculate Total Project Quality based on Quality Domain for a project 
%structure
%----------------
%Output:
%TPQ: Total Project Quality (scalar)
%---------------- 
%Inputs:
%DSM: Upper triangular binary matrix of logic domain (a project structure
%of a PEM matrix)
%PEM: Upper triangular matrix of logic domain
%QD:  N by w matrix of quality domain
%q:   N by 1 vector of quality parameters
%---------------- 
%Usage:
%TPQ=tpqfastqd(DSM,PEM,QD,q)
%---------------- 
%Example: Calculate TPQ for a (best) project structure 
%PEM=triu(rand(10)*.5+.5); %Generate PEM
%DSM=round(PEM); %Genarate DSM
%QD=rand(10,2); Generate QD domain
%q=QD(:,1); %Specify q vector
%TPQ=tpqfastqd(DSM,PEM,QD,q) %Calculate TPQ
%---------------- 
%Prepositions and Requirements:
%1.)The diagonal values of PEM must be in interval [0,1]
%2.)The diagonal values of DSM must be binary values
%3.)The elements of column vector q and QD matrix must be in interval [0,1]

function TPQ=tpqfastqd(DSM,PEM,QD,q)
pem=diag(PEM);
dsm=diag(DSM);
TPQ=0;
if sum(max(QD(pem>0,:),[],2))>0
    TPQ=sum(q(dsm>0))/sum(max(QD(pem>0,:),[],2));
end
