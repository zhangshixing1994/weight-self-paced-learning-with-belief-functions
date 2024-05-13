function [Xtrain,ytrain,Xtest,ytest] = load_data


data=load('data.txt');
X=data(:,2:5);
label=ones(50,1)*[1:3];
y=label(:);
nClass=3;

[X,y,nClass] = NewDataSet(1);

X=data_standard(X);
y_=onehot(y);

[idx_train,idx_test] =train_test_split2(y,0.5);

Xtrain=X(idx_train,:);ytrain=y(idx_train);
Xtest=X(idx_test,:);ytest=y(idx_test);ytest_=y_(idx_test,:);

[ntrain,dim]=size(Xtrain);
[ntest,~]=size(Xtest);
Xtrain = [ones(ntrain,1), Xtrain]; 
Xtest=[ones(ntest,1),Xtest];

end

