%Author: Zsolt T. Kosztyán Ph.D habil., University of Pannonia,
%Faculty of Economics, Department of Quantitative Methods
%----------------
%Calculate Total Project Quality for a project structure
%----------------
%Output:
%TPQ: Total Project Quality (scalar)
%---------------- 
%Inputs:
%DSM: Upper triangular binary matrix of logic domain (a project structure
%of a PEM matrix)
%PEM: Upper triangular matrix of logic domain
%q:   N by 1 vector of quality parameters
%---------------- 
%Usage:
%TPQ=tpqfast(DSM,PEM,q)
%---------------- 
%Example: Calculate TPQ for a (best) project structure 
%PEM=triu(rand(10)*.5+.5); %Generate PEM
%DSM=round(PEM); %Genarate DSM
%q=rand(10,1); %Generate q vector (quality parameters)
%TPQ=tpqfast(DSM,PEM,q) %Calculate TPQ
%---------------- 
%Prepositions and Requirements:
%1.)The diagonal values of PEM must be in interval [0,1]
%2.)The diagonal values of DSM must be binary values
%3.)The values of column vector q must be in interval [0,1]

function TPQ=tpqfast(DSM,PEM,q)
TPQ=0; %Total Project Quality
TPS=0; %Total Project Score (additive scores are assumed)
for i=1:size(DSM,1)
    if DSM(i,i)>0
        TPS=TPS+PEM(i,i); %TPS=S(PEMii)
        if TPQ==0
            TPQ=1;
        end
        TPQ=TPQ*power(q(i),PEM(i,i)); 
    end
end
if TPS>0
    TPQ=maxscore_PEM(DSM,PEM,ones(size(PEM,2))-PEM)*power(TPQ,1/TPS); 
    %TPQ=SCOPE*TTi(qi^PEMii)^(1/TPS)
end