
function [net] = train_model(net,Xtrain,ytrain,weight)


input_layer_size=net.input_layer_size;    
hidden_layer_size=net.hidden_layer_size;
max_iter=net.max_iter;
alpha=net.alpha;
[ntrain,nClass]=size(ytrain);



J_history=zeros(max_iter,1);

w1=net.w1;
w2=net.w2;

iter = 1;

while iter<max_iter
% Forward Propagation

z2 = Xtrain*w1'; %
a2 = sigmoid(z2); %

a2 = [ones(ntrain,1) a2];
z3 = a2*w2';
od = sigmoid(z3); 

%logistic regression
J = 1/ntrain * sum(sum(-ytrain .* log(od) - (1 - ytrain) .* log(1 - od)));
J_history(iter)=J; 

% Backward Propagation

%delta_3 = od - ytrain;
delta_3 = (od - ytrain).*repmat(weight,1,nClass);

delta_2 = ((delta_3 * w2(:,2:end)) .* sig_grad(z2));

w1_grad = 1/ntrain* delta_2' * Xtrain;
w2_grad = 1/ntrain * delta_3' * a2;

%Weight Update
grad=[ w1_grad(:);w2_grad(:)];

w=[w1(:);w2(:)];

w=w-alpha*grad;

w1=reshape(w(1:hidden_layer_size * (input_layer_size + 1)), hidden_layer_size, (input_layer_size + 1));
w2=reshape(w(1 + (hidden_layer_size * (input_layer_size + 1)):end),nClass, (hidden_layer_size + 1));

iter=iter+1;

end
net.w1=w1;
net.w2=w2;
end