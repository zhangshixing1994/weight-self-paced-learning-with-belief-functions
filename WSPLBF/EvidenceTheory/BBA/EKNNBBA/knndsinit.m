function [gamm,alpha] = knndsinit(x,S);

% Copyright Thierry Denoeux 
% February 14, 2001
%
% KNNDSINIT: [gamma,alpha] = knndsinit(x,S)
%	
% 	This function initialises parameter gamma and alpha of the BPA 
%	for every class
% 	gamma (q) is the inverse of the mean distance
%	between two vectors of class q
%
%	gamma: vector (M,1) of parameters of the BPA.
%	 M is the number of classes.
%	alpha (1,1): also a parameter of the BPA.
% 	x: matrix (N,d) of the training set patterns
% 	S: vector (N,1) of the training set labels
%
%	See also: KNNDSFIT,KNNDSVAL
%
% References:
% 
% T. Denoeux. A k-nearest neighbor classification rule based on 
%  Dempster-Shafer theory. IEEE Transactions on Systems, Man
%  and Cybernetics, 25(05):804-813, 1995.
%
% L. M. Zouhal and T. Denoeux. An evidence-theoretic k-NN rule with 
% parameter optimization. IEEE Transactions on Systems, Man and 
% Cybernetics - Part C, 28(2):263-271,1998.


[Napp,nent]=size(x);
M=max(S);


for i=1:M,
  ii=find(S==i);Nii=length(ii);
  D=zeros(Nii,Nii);
  for j=1:Nii
   D(j,:)=sqrt(sum(((ones(Nii,1)*x(ii(j),:))-x(ii,:))'.^2)')';
   end;
	%moyenne de ces Nii(Nii-1) distances
Dm(i) = sum(sum(D))/(Nii*Nii - Nii);
end;

gamm = ones(1,M) ./ Dm;
gamm=gamm';
alpha=.95;

