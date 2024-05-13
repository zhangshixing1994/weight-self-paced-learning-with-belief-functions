

function [output] = WSPLBF(Xtrain,ytrain,Xtest,ytest,net,AM,xi,weight)
nClass=length(unique(ytest));
ntrain=length(ytrain);
ytrain_=onehot(ytrain);
AM=AM/max(AM);





V = zeros(ntrain,1);

for j=1:6
    
  V(:) = 0;
  output= MyBP_test(net,Xtrain,ytrain);
  E=output.E/max(output.E);




  for k=1:nClass
      Loss=(1-xi(j))*E+xi(j)*AM;
      idx = find(ytrain==k);
      m_k = length(idx);

      sortLoss_k = sort((1-xi(j))*E(idx)+xi(j)*AM(idx));
      lambda_SPL_k = sortLoss_k(round(m_k*(0.1+0.15*j)));
      tmp = Loss<=lambda_SPL_k;
      tmp(ytrain~=k) = 0;
      V(tmp) = 1;
            
  end



net = train_model(net,Xtrain(V==1,:),ytrain_(V==1,:),weight(V==1));




output = inference(net,Xtest,ytest);


  
        
        
           
end   


output.net=net;


end