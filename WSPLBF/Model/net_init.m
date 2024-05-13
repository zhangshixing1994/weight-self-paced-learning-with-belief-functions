function [net] = net_init(dim,nClass,max_iter,alpha)


%ÍøÂç²ÎÊıÉèÖÃ
net.input_layer_size  = dim;   
net.hidden_layer_size = 2*dim;   
net.max_iter=max_iter;
net.alpha=alpha;
w1 = random(net.input_layer_size, net.hidden_layer_size);
w2 = random(net.hidden_layer_size, nClass);
net.w1=w1;
net.w2=w2;
end

