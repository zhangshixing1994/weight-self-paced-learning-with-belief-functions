function [AM] = calc_AM(Xtrain,ytrain,Xtest,k)

[mass] = EKNNBBA(Xtrain,ytrain,Xtest,k);

% [n,m]=size(mass);
% BetP=mass(:,1:m-1);
% for i=1:n
%     BetP(i,:)=BetP(i,:)/(1-mass(i,m));
%     
% end
% AM=-sum(BetP.*log(BetP)/log(2),2);

[AM,BetP] = compute_AM(mass);


end

