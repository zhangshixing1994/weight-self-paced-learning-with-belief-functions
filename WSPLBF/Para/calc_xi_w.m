function [xi,w] = calc_xi_w(Xtrain,ytrain,Xtest,ytest,net,AM,w)


psi=1;
xi=[1,0.8,0.6,0.4,0.2,0];% training rounds is set as 6 in our datasets;if you use large dataset,xi needs adjust.
[bestParam] = search_para(Xtrain,ytrain,Xtest,ytest,net,AM,xi,w);
gama= bestParam.gama;
psi=bestParam.psi;

xi=xi.^psi;
w=w.^gama;
end

