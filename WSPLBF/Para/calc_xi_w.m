function [xi] = calc_xi(Xtrain,ytrain,Xtest,ytest,net,AM)


psi=1;
xi=[1,0.8,0.6,0.4,0.2,0];
[bestParam] = search_para(Xtrain,ytrain,Xtest,ytest,net,AM,xi);
gama= bestParam.gama;
psi=bestParam.psi;

xi=xi.^psi;
end

