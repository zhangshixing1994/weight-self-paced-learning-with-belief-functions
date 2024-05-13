

function [mass] = EKNNBBA_train(Xtrain,ytrain,k)


[gamm,alpha] = knndsinit(Xtrain,ytrain);

ntrain=length(ytrain);






[m,~] = knndsval1(Xtrain,ytrain,k,gamm,alpha,0,Xtrain);


[p,q]=size(m);
mass=zeros(p,2^(q-1));
for i=1:q-1
    mass(:,2^(i-1)+1)=m(:,i);
end
mass(:,2^(q-1))=m(:,q);
end