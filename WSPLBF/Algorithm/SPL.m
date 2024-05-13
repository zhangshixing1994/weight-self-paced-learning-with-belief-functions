function [output] = SPL(Xtrain,ytrain,Xtest,ytest,net)
nClass=length(unique(ytest));
ntrain=length(ytrain);
ytrain_=onehot(ytrain);

weight=ones(length(ytrain),1);
V0=[1:1:ntrain];
V = zeros(ntrain,1);
Acc=[];
for j=1:6
    
  V(:) = 0;


output= inference(net,Xtrain,ytrain);


  
  for k=1:nClass
            idx = find(ytrain==k);
            m_k = length(idx);
            sortLoss_k = sort(output.E(idx));
            lambda_SPL_k = sortLoss_k(round(m_k*(0.1+0.15*j)));
            tmp = output.E<=lambda_SPL_k;
            tmp(ytrain~=k) = 0;
            V(tmp) = 1;
  end


net = train_model(net,Xtrain(V==1,:),ytrain_(V==1,:),weight(V==1));


output = inference(net,Xtest,ytest);

end   


output.net=net;


end