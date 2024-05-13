function [gamm,alpha,err] = knndsfit(x,S,K,gamm,version,param);

% Copyright Thierry Denoeux 
% February 14, 2001
%
% KNNDSFIT:[gamm,alpha,err] = knndsfit(x,S,K,gamm,version,param)
%
%	provides the training error rate and fits the parameters 
%	of the K-nearest neighbour classification rule based on 
%	Dempster-Shafer theory
%
% 	Inputs:
%
%	K : number of neighbours
% 	x : matrix (napp,d) of the training set 
% 	S : vector (napp,1) of corresponding labels
% 	version: 0 : without optimising the parameters of the BPA
%	   	 1 : version optimising the parameters by an 
%	algorithm of gradient descent 
% 	gamm: vector (M,1) of parameters of the BPA
%	param : optional vector of parameters for both versions
%	(usual values have been fixed by default)  
%	version 0 : alpha, parameter of the BPA  (0.95)
%	version 1 : (1,4) vector containing:
% 		* alpha (0.95)
%		* maximum number of iterations in the
%		 optimisation loop (100) 
%               * initial step of gradient variation (0.1)
%%		* minimum gain in the optimisation loop (1e-6)
%
% 	Outputs:
%		
% 	err : training error rate
%	gamm: vector (M,1) of optimised parameters of the BPA
%
%	The method can be divided in 3 steps:
%
%	* for every x of the training set ,computing the K-nearest
%	 neighbours
%	* optionnally, optimising the parameters of the BPA
%	* computing the error rate
%
%	See also: KNNDSINIT,KNNDSVAL
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


if nargin < 3
	disp('this function needs at least 3 arguments')
end;

if nargin== 3
	version=0;
	[gamm,alpha]=knndsinit(x,S);
	param(1)=alpha;
end;
	
if nargin== 4
	version=0;
	alpha=.95;
	param(1)=alpha;	
end;

if nargin == 5
	if version ==1
		param(1)=.95;
		param(2)=100;
		param(3)=.1;
		param(4)=1e-6;
		alpha=param(1);
	else
		alpha=.95;
	end;
end;


if nargin==6
  	alpha=param(1);
end;

if nargin >=7 
	disp('too many input arguments')
end;

if nargout >=4 
	disp('too many output arguments')
end;

[Napp,nent]=size(x);

ncl=max(S);
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul des K ppv dans l'ens. d'app.

ds=[];is=[];
for i = 1:Napp,
	dist=sum(((ones(Napp,1)*x(i,:))-x)'.^2)';
	[dss,iss]=sort(dist);
	ds = [ds dss(2:K+1)];
	is = [is iss(2:K+1)];
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determination des parametres de la BPA

% choix des parametres de la methode d'optimisation ou option par defaut

if version==1,
	gamm=gamm';
	gamm = optimds(x,S,gamm,ds,is,K,param);
 end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul de l'erreur de test et d'apprent.

 
err = classds(alpha,gamm,ds,is,S,K);

function  gamm = optimds(x,S,gamm,ds,is,K,param)

      
% OPTIMDS: Algorithm optimising the parameters of the Basic 
%	Probability Assignment by gradient descent of the 
%	classification error
%	gamm = optimds(x,S,ds,is,K,param)
%
%	x: matrix (N,d) of the training set patterns
% 	S: vector (N,1) of the training set labels
% 	ds,is: matrices (K,N) containing distances and indices of 
%	the K-nearest neighbours of N vectors of the training set
%	K: number of neighbours
%	param: vector (1,4) containing:
%		(usual values are given)
% 		* the parameter alpha of the BPA (0.95)
%		* maximum number of iterations in the
%		 optimisation loop (100) 
%               * initial step of gradient variation (0.1)
%		* minimum gain in the optimisation loop (1e-6)
%
% 	gamm: vector (M,1) of optimised parameters 
%		(M is the number of classes)
%
%	See also: KNNDS, CLASSDS, GRADIENTDS, INITDS

M = max(S);
N = length(S);
Ident = eye(M);
T = Ident(:,S);

alpha =param(1);maxiter=param(2);
eta =param(3);a =1.2;b=.8;c = .5;
mi = 1e-4;mx = 1e6;gain_min=param(4);

pas = eta * ones(M,1);
it = 0; gain = 1;


P=gamm';


[Errcou1,Dp] = gradientds(S,ds,is,P,alpha,K);
Pp = P;
Errcop = Errcou1 + 1;

while gain >= gain_min & it <= maxiter,
	it = it + 1;
	[Errco,D] = gradientds(S,ds,is,P,alpha,K);
	if isnan(gain) | isinf(gain), gain = 1; end;
	if Errco > Errcop,
		P = Pp;
		D = Dp;
		pas = pas .* c;
		P = P - pas .* D;
	else
		gain = .9*gain + .1*abs(Errcop - Errco);
		Pp = P;
		test =  (D .* Dp) >= 0; 
		pas = ((test * a) + ((~test) * b)) .* pas;
		pas = (pas <= mx) .* pas + (pas > mx) * mx;
		pas = (pas >= mi) .* pas + (pas < mi) * mi;
		Dp = D;
		P = P - pas .* D;
		Errcop = Errco;
	end;
	
	gamm= P;
 end;



function [ERR,D] = gradientds(S,ds,is,P,alpha,K);

% GRADIENTDS:
%	[ERR,D] = gradientds(S,ds,is,P,alpha,K)
%	calculates the error and gradient vector 
%	in the algorithm optimising the gamma parameters
%	of the BPA in the KNNDS method of classification	
%	 
% 	S: vector (N,1) of labels of the training set
% 	ds,is: matrices (K,N) containing distances and indices of 
%	the K-nearest neighbours of N vectors of the training set
%	K: number of neighbours
%	P : vector (1,M) of squared parameters gamma of the BPA
%	alpha: parameter of the BPA
%	(M is the number of classes)

% 	ERR:  classification cost
% 	D: gradient vector (1,M)
%
%	See also: OPTIMDS, KNNDS

N=length(S);
M = max(S);
Ident = eye(M);
T = Ident(:,S);


gama = P(1:M); 

Dgama = zeros(M,1);
Ds = zeros(1,N);
Dsgama = zeros(M,N);

mk = [zeros(M,N);ones(1,N)];
mm = mk;
s = zeros(K,N);
ss = zeros(K,N);
lamda = 1/M;

for k = 1:K,
 	Is = is(k,:);
	Is = S(Is);
	G = zeros(M,N);
	Tk = zeros(M,N);
   for j = 1:M,
      pos = find(Is==j);
	if length(pos) ~= 0
	   Tk(j,pos) = ones(1,length(pos));
	end;
 	end;
	G = gama.^2 * ones(1,N).*Tk;
	gam = max(G);
   ss(k,:) = exp(-gam .*ds(k,:));
	s(k,:) = alpha*ss(k,:);
 	m = [Tk.*(ones(M,1)*s(k,:));1-s(k,:)];
	mk = [mk(1:M,:).*(m(1:M,:)+(ones(M,1)*m(M+1,:)))+...
	m(1:M,:).*(ones(M,1)*mk(M+1,:));(mk(M+1,:).*m(M+1,:))];
end;

%             normalisation

Kn = sum(mk);
mkn = mk./(ones(M+1,1)*Kn);
Q = mkn(1:M,:)+lamda*ones(M,1)*mkn(M+1,:) - T;

ERR = 0.5*mean(sum(Q.^2));

%		gradient

Dsm = zeros(M+1,N);
v = zeros(M,N);

for k = 1:K
	Is = is(k,:);
	Is = S(Is);
   Tk = zeros(M,N);
   for j = 1:M,
      pos = find(Is==j);
	if length(pos) ~= 0
		Tk(j,pos) = ones(1,length(pos));
	end;
 	end;
	m = [Tk.*(ones(M,1)*s(k,:));1-s(k,:)];
	if length(find(m(M+1,:) == 0)) > 1 ,
	   Dsgama = 1e10*ones(M,N);
	   k = K+1;
	else
	   mm(M+1,:) = mk(M+1,:)./m(M+1,:);
	   H = ones(M,1)*mm(M+1,:);
	   mm(1:M,:) = (mk(1:M,:) - H.*m(1:M,:))./(m(1:M,:)+ones(M,1)*m(M+1,:));
	   v = zeros(M,N);
	   v(1:M,:) = (mm(1:M,:) + ones(M,1)*mm(M+1,:)).*Tk - mm(1:M,:);
	   DsK = sum(v) - mm(M+1,:);
	   Dsm(1:M,:) = ((ones(M,1)*Kn).*v - mk(1:M,:).*(ones(M,1)*DsK))./(ones(M,1)*Kn.^2);
	   Dsm(M+1,:) = (-Kn.*mm(M+1,:) - mk(M+1,:).*DsK)./(Kn.^2);	
	   Ds = sum(Q.*(Dsm(1:M,:)+lamda*ones(M,1)*Dsm(M+1,:)));
	   Dsgama = Dsgama -2*gama*(Ds.*ds(k,:).*s(k,:)).*Tk;   
 	end;
end;
Dgama = ((mean((Dsgama)'))'); 
D = Dgama;

function [err,m,L] = classds(alpha,gamm,ds,is,S,K);

% CLASSDS:  [err,m,L] = classds(alpha,gamm,ds,is,S,K)
%
%	This program provides the BPA, the labels and the
%	error rate  in a K-nearest neighbour rule based on 
%	Dempster-Shafer theory.  
% 		
% 	alpha (1,1), gamm (1,M): parameters of the BPA
%	(M is the number of classes)
% 	ds, is : matrices (K,N) containing respectively 
%	distances et indices of the K - nearest neighbours 
%	of N vectors belonging to the test set.
%	S: vector (N,1) of labels of the test set. 
%	S: vector (napp,1) of labels of the training set.
%	K: number of neighbours
%
% 	L: vector (N,1) of final labels
% 	err: test error rates
% 	m: matrix (M+1,N) of the final BPA 

N=length(S);
ncl=max(S);
m = [zeros(ncl,N);ones(1,N)]; % element neutre de la combinaison de Dempster

cppv=zeros(N,1);
	
for i=1:N,
   for j=1:K,
     m1 = zeros(ncl+1,1);
     m1(S(is(j,i))) = alpha*exp(-gamm(S(is(j,i))).^2*ds(j,i));
     m1(ncl+1) = 1 - m1(S(is(j,i))); 
     m(1:ncl,i) = m1(1:ncl).*m(1:ncl,i) + m1(1:ncl)*m(ncl+1,i) + m(1:ncl,i)*m1(ncl+1);
     m(ncl+1,i) = m1(ncl+1) * m(ncl+1,i);    
  end;
end;
m=m./(ones(ncl+1,1)*sum(m));
[temp,L]=max(m(1:ncl,:));
L=L';
err=length(find(L ~= S))/N;


