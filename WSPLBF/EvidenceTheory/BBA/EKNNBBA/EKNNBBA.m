function [mass] = EKNNBBA(train_data,train_label,test_data,k)


[gamm,alpha] = knndsinit(train_data,train_label);
[m,L] = knndsval(train_data,train_label,k,gamm,alpha,0,test_data);

[p,q]=size(m);
mass=zeros(p,2^(q-1));
for i=1:q-1
    mass(:,2^(i-1)+1)=m(:,i);
end
mass(:,2^(q-1))=m(:,q);
end

