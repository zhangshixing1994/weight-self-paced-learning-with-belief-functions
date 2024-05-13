function [output] = inference(net,Xtest,ytest)

w1=net.w1;
w2=net.w2;
ntest=length(ytest);

z2 = Xtest * w1';
a2 = [ones(size(sigmoid(z2), 1), 1) sigmoid(z2)];

% H_theta(x)
z3 = a2 * w2';
p = sigmoid(z3);

[~, ypred] = max(p, [], 2);
acc=length(find(ypred==ytest))/ntest;

ytest_=onehot(ytest);
E=sum(-ytest_ .* log(p) - (1 - ytest_) .* log(1 - p),2);

output.p=p;
output.ypred=ypred;
output.acc=acc;
output.E=E;
end