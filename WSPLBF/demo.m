

clc; clear all; close all;
addpath(genpath(pwd))



[Xtrain,ytrain,Xtest,ytest] = load_data;
nClass=length(unique(ytrain));
[~,dim]=size(Xtrain);dim=dim-1;

net = Init_model(dim,nClass);

output = SPL(Xtrain,ytrain,Xtest,ytest,net);

[AM,w] = calc_evidinfo(Xtrain,ytrain);
[xi,w] = calc_xi_w(Xtrain,ytrain,Xtest,ytest,net,AM,w);
output2 = WSPLBF(Xtrain,ytrain,Xtest,ytest,net,AM,xi,w);