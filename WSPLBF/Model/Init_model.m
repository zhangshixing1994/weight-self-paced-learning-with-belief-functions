function [net] = Init_model(dim,nClass)

max_iter=1000;
alpha=0.2;
net = net_init(dim,nClass,max_iter,alpha);

end

