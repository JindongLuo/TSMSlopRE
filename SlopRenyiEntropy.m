function [SlopEn,npdf] = SlopRenyiEntropy(x, m,delta,gama,q)
%   Returns the slope entropy estimate of the data sequence 

%      m : Embedding Dimension, an integer > 1. 
%      delta : threshold
%      gama : threshold
%      gama > delta 
if gama<delta
    disp('gama should be greater than delta ')
end
N=length(x);
for j=1:N-m+1
    for i=j+1:j+m-1        
        if -delta<=x(i)-x(i-1) && x(i)-x(i-1)<=delta
            Symbx_temp(i-j)=0;
        elseif delta<x(i)-x(i-1) && x(i)-x(i-1)<=gama
            Symbx_temp(i-j)=1;
        elseif gama<x(i)-x(i-1)
            Symbx_temp(i-j)=2;
        elseif -gama<=x(i)-x(i-1) && x(i)-x(i-1)<-delta
            Symbx_temp(i-j)=-1;
        elseif x(i)-x(i-1)<-gama
            Symbx_temp(i-j)=-2;
        end
    end
    Symbx(j,:)=Symbx_temp;    % symbolic x
end
% all symbol patterns
all_patterns=[-2:2]';
for f=2:m-1
    temp=all_patterns;
    all_patterns=[];
    j=1;
    for w=1:length(temp)
        [a,b]=size(temp);
        all_patterns(j:j+a-1,:)=[temp,temp(w)*ones(a,1)];
        j=j+a;
    end
end
[all_patterns,~,~]=unique(all_patterns,'rows');
% statistics
npdf(1:size(all_patterns,1))=0;
 for i=1:N-(m-1)
     
     for jj=1:size(all_patterns,1)
         if (abs(all_patterns(jj,:)-Symbx(i,:)))==0
             npdf(jj) = npdf(jj) + 1 ;
         end
     end
 end
% calculate
npdf= npdf./(N -(m-1));
p=npdf(npdf~=0); % Remove cases where probability is 0
% -------------------------------------------------------------------------
SlopEn =log(sum(p.^q))/(1-q); %%RenyiìØ
end