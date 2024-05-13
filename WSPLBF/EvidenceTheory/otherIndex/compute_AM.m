function [AM,BetP] = compute_AM(Mass,ytrain)


%����������BetP
n=size(Mass,1);
BetP=[];
for i=1:n
betp=mtobetp(Mass(i,:));
BetP=[BetP;betp];
end
%��������AM
AM=-sum(BetP.*log(BetP)/log(2),2);

[~,ypre]=max(BetP');

idx=find(ypre'~=ytrain);
AM(idx)=max(AM);
end

