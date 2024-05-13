function [AM,BetP] = compute_AM(Mass,ytrain)


%计算样本的BetP
n=size(Mass,1);
BetP=[];
for i=1:n
betp=mtobetp(Mass(i,:));
BetP=[BetP;betp];
end
%计算样本AM
AM=-sum(BetP.*log(BetP)/log(2),2);

[~,ypre]=max(BetP');

idx=find(ypre'~=ytrain);
AM(idx)=max(AM);
end

