function [Mass] = ECMBBA(train_data,train_label,test_data,beta,alpha,delta)
[gplus,c] = Class_center(train_data,train_label);
k=length(unique(train_label));
[Mass] = calculate_mass(test_data,gplus,c,k,beta,alpha,delta);
end

