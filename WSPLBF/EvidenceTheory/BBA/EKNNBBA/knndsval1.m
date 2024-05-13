

function [m,L] = knndsval1(xtrain,ytrain,K,gamm,alpha,loo,xtest);


[ntrain,nent]=size(xtrain);
M=max(ytrain);

if loo,
   xtest=xtrain;
end;

[ntest,nent]=size(xtest);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of the K-nearest neighbours in the training set

  dst=[];ist=[];
  for i = 1:ntest,
	dist=sum(((ones(ntrain,1)*xtest(i,:))-xtrain)'.^2)';
	[dss,iss]=sort(dist);
   dst = [dst dss(2+loo:K+loo)];
   ist = [ist iss(2+loo:K+loo)];
  end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of the BPA

m = classdstst(alpha,gamm,xtest,dst,ist,ytrain,K); 
[temp,L]=max(m(:,1:M)');
L=L';



%%%%%%%%%%%% dependent programs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function m = classdstst(alpha,gamm,xtest,dist,is,ytrain,K);

% CLASSDSTST:  m = classds(alpha,gamm,ds,is,Sapp,K)
%
%	This program provides the BPA, the labels and the
%	error rate  in a K-nearest neighbour rule based on 
%	Dempster-Shafer theory.  
% 		
% 	alpha (1,1), gamm (M,1): parameters of the BPA
%	(M is the number of classes)
% 	ds, is : matrices (K,N) containing respectively 
%	distances et indices of the K - nearest neighbours 
%	of N vectors belonging to the test set. 
%	Sapp: vector (napp,1) of labels of the training set.
%	K: number of neighbours
%
% 	m: matrix (M+1,N) of the final BPA 

N= max(size(xtest));
M=max(ytrain);
m = [zeros(M,N);ones(1,N)]; 
cppv=zeros(N,1);
	
for i=1:N,
   for j=1:K-1
     m1 = zeros(M+1,1);
     m1(ytrain(is(j,i))) = alpha*exp(-gamm(ytrain(is(j,i))).^2*dist(j,i));
     m1(M+1) = 1 - m1(ytrain(is(j,i)));
     m(1:M,i) = m1(1:M).*m(1:M,i) + m1(1:M)*m(M+1,i) + m(1:M,i)*m1(M+1);
     m(M+1,i) = m1(M+1) * m(M+1,i);    
  end;
end;
m=m./(ones(M+1,1)*sum(m));
m=m';