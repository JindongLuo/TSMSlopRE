   function [TSEn] = TSMSlopRE(X,Kmax,dim,delta,gama,q)

%{
Time-Shift Multiscale Slop Entropy (TSMSlopEn): compute average entropy values for Kmax
Script for computing predictability of irregular time series at multiple intervals of time series.
Notion is based on Higuchi fractal dimension

Functions called: TSE_beta.m

INPUT:
    S: original time series
    Kmax: maximum time intervals
    dim: embedding dimension or template length
    r: tolerance value to be multipled with SD
 
    
OUTPUT:
    TSEn: vector of average K values of entropy

%}

TSEn=zeros(1,Kmax);

for k=1:Kmax
    [En] = TSE_beta(X,k,dim,delta,gama,q);
    TSEn(k) = mean(En); % mean En
end
end

function [En] = TSE_beta(S,K,dim,delta,gama,q)

%{
Time-Shift Entropy (TSE)
Script for computing predictability of irregular time series at multiple intervals of time series.
Notion is based on Higuchi fractal dimension,

INPUT:
    S: original time series
    K: time interval
    dim: embedding dimension or template length
    r: tolerance value to be multipled with SD
    
OUTPUT:
    En: vector of K (beta's) values of entropy

%}

% Composing of the sub-series:
N = length(S);
X{K}=[];

for beta = 1:K
    for m = 1:beta
        X{beta}=[];
        limit = floor((N-m)/beta);
        for i = m:K:(m + (limit*beta))
            X{beta} = [X{beta},S(i)];
        end
    end
end

% compute Slop Entropy for each sub-time series
for i=1:K
    [En(i),~] = SlopRenyiEntropy(X{i},dim,delta,gama,q);% Slop Entropy
end
end


