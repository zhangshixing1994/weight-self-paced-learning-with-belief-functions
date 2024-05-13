function [m,g,F,pl,BetP,J,N] = ECM(x,K,version,alpha,beta,delta,init);
% Evidential C-Means
%   [m,g,F,pl,BetP,J,N] = ECM(x,K,version,alpha,beta,delta,init);
% INPUTS
%   x: input matrix nxd
%   K: number of desired clusters
%   Optional : 
%       version : 0 (default):  2^K focal elements, 1: focal elements of
%       size less or equal to 2 except Omega 
%       alpha : exponent for the cardinality (default 1) 
%       beta : exponent for m (defaut 2) 
%       delta : distance to the empty set (default 10)
%       init : initialization method 0 (default) random, 1 : by fuzzy c-means  
%       
% OUTPUTS
%   m: credal partition = maxtrix of size nxf; the corresponding focal elements are given in F; 
%   g: matrix of size Kxd of the centers of the clusters
%   F: matrix (f,K) of focal elements
%       F(i,j)=1 if omega_j belongs to focal element i
%              0 otherwise
%   pl: plausibilities of the clusters
%   BetP: pignistic probabilities of the clusters
%   J: objective function
%   N: non specifity index (to be minimized)
%
% Reference:
% M.-H. Masson and T. Denoeux. ECM: An evidential version of the fuzzy c-means algorithm. 
% Pattern Recognition, Vol. 41, Issue 4, pages 1384– 1397, 2008.
%
% June 8, 2011.


% --------------- controls and default values --------------------
[n,d]=size(x);
if nargin<2
    error('ECM needs two arguments');
end

if nargin==2
version=0;
beta=2;alpha = 1;
delta = 10 ;  init=0;
end

if nargin==3
beta=2;alpha = 1;
delta = 10 ;  init=0;
end

if nargin==4
beta=2;
delta = 10 ;  init=0;
end

if nargin==5
delta = 10 ; init=0;
end

if nargin==6
 init=0;
end

if isempty(alpha) alpha=1;end
if isempty(beta) beta=2;end
if isempty(delta) delta=10;end
if isempty(version) version=0;end
if isempty(init) init=0;end
fprintf('-------------------------------------------\n');
fprintf('Evidential c-means\n');
fprintf('-------------------------------------------\n');
fprintf('Number of objects = %5i\n',n);
fprintf('alpha= %5.2f\n',alpha);
fprintf('beta = %5.2f\n',beta);
fprintf('delta = %5.2f\n',delta);

%--------------------------------------------------------------
delta2=delta^2;
%--------------- construction of the focal set matrix  ---------
ii=1:2^K;
F=zeros(length(ii),K);
for i=1:K,
      F(:,i)=bitget(ii'-1,i);
end;

if version==0 nbFoc=2^K;end

if version==1  % limitation of focal sets to cardinality <=2 + Omega

    if K>3
    truc = sum(F,2);
    ind = find(truc>2);
    ind(end)=[];
    F(ind,:)=[];
    end

    FF = F ; 

    if K~=2
        nbFoc = K + K*(K-1)/2 + 2 ; % with empty set
    else
        nbFoc=4;
    end;
end;

c = sum(F(2:end,:),2)';  % cardinality of focal sets
if version==1 fprintf('The cardinality of the focal elements limited\n');end
fprintf('Number of focal elements= %2i\n',nbFoc);
    
%---------------------- initialisations --------------------------------------
if init ==0
     fprintf('Random initialization of the centers\n');
g = randn(K,d);
else
    fprintf('Initialization of the centers using FCM\n');
    g = FCM(x,K);
end

pl=[];
BetP=[];
histJ=[];
%------------------------ iterations--------------------------------
fprintf('-----------------------------------------\n');
fprintf('Optimization\n');
pasfini=1;Jold = inf; 
while pasfini

gplus=[];
for i=2:nbFoc
    fi = F(i,:);
    truc = repmat(fi',1,d);
    gplus = [gplus;sum(g.*truc)./sum(truc)];
end

    % calculation of distances to centers 
    D=[];
    for j=1:nbFoc-1
        A = (x - ones(n,1)*gplus(j,:));
        B = diag(A*A');
        D = [D B];
    end
       
    % Calculation of masses
    
    m = zeros(n,nbFoc-1);
    for i=1:n
        vect0 = D(i,:);
        for j=1:nbFoc-1
            vect1 = ((D(i,j)*ones(1,nbFoc-1))./vect0).^(1/(beta-1));
            vect2 =  ((c(j)^(alpha/(beta-1)))*ones(1,nbFoc-1))./(c.^(alpha/(beta-1)));
            vect3 = vect1.*vect2 ;
            m(i,j)=1/(sum(vect3)+(c(j)^alpha*D(i,j)/delta2)^(1/(beta-1)));
        end
    end
   
    % Calculation of centers
    
    A = zeros(K,K);
    B = zeros(K,1);
   
    for k=1:K
        for l=1: K
        truc = zeros(1,K);
        truc(1,k)=1;truc(1,l)=1;
        t = repmat(truc,nbFoc,1) ; 
        indices = find(sum((F-t)-abs(F-t),2)==0) ;   % indices of all Aj containing wk and wl
        indices = indices - 1 ;
        for jj = 1:length(indices)
            j = indices(jj);
            mj = m(:,j).^beta;
            A(l,k)=A(l,k)+sum(mj)*c(j)^(alpha-2);
        end
        end
    end
    
    B=[];
    for l=1:K
        truc = zeros(1,K);
        truc(1,l)=1;
        t = repmat(truc,nbFoc,1) ; 
        indices = find(sum((F-t)-abs(F-t),2)==0) ;   % indices of all Aj containing wl
        indices = indices - 1 ;
        mi = repmat((c(indices).^(alpha-1)),n,1).*(m(:,indices).^beta);
        s = sum(mi,2);
        mats = repmat(s,1,d);
        xim = x.*mats ;
        blq = sum(xim);
        B=[B;blq];
    end
    

    g=A\B;
    
    mvide = ones(n,1)-sum(m,2) ; 
    J = sum(sum((m.^beta).*D.*repmat(c.^alpha,n,1)))+ delta2*sum(mvide.^beta);
    fprintf('Objective function = %6.3f\n',J);
    pasfini =abs(J-Jold)>1e-3;
    Jold = J ; 
    histJ=[histJ;J];
end
%--------------------------- end of iterations ----------------------------
m = [ones(n,1)-sum(m,2) m]; 

if version==1
    if K>3
        mm=zeros(n,2^K);
        ii=1:2^K;
        F=zeros(length(ii),K);
        for i=1:K,
            F(:,i)=bitget(ii'-1,i);
        end;
        truc = sum(F,2);
        ind = find(truc<=2);ind=[ind;2^K];
        for j=1:length(ind)
            mm(:,ind(j))=m(:,j);
        end
    else
        mm=m;
    end
    P=[];
    for i=1:n
        pp = mtopl(mm(i,:));
        P=[P;pp];
        bet = betp(mm(i,:));
        BetP = [BetP;bet];
    end
    F=FF;
else
P=[];
for i=1:n
    pp = mtopl(m(i,:));
    P=[P;pp];
    bet = betp(m(i,:));
    BetP = [BetP;bet];
end
end

truc = sum(F,2);
singletons = find(truc==1);
pl = P(:,singletons);

% ----------- validity index ---------------------------
truc = truc';
mat = repmat(truc,n,1);
mat(:,1)=K*ones(n,1);
N = sum(sum(log(mat).*m))/log(K)/n;

fprintf('End of optimization\n');
fprintf('-------------------------------------------\n');


function [out] = betp(in) 
% computing BetP on ? from the m vector (in) 
% out = BetP vector: order a,b,c,... 
% beware: not optimalize, so can be slow for >10 atoms 
 
lm = length(in); 
natoms = round(log2(lm)); 		 
if 2^natoms == lm  
	if in(1) == 1 
		out = 	ones(1,natoms)./natoms; 
	else 
		betp = zeros(1,natoms); 
		for i = 2:lm 
			x = bitget(i-1,1:natoms); % x contains the 1 and 0 of i-1, for a,b,c... 
			betp = betp + in(i)/sum(x).*x; 
		end 
		out = betp./(1-in(1)); 
	end 
else 
	'ACCIDENT in betp: length of input vector not OK: should be a power of 2' 
end 


function [out] = mtopl(in) 
% computing FMT from m to pl 
% in = m vector 
% out = pl vector 
 
in = mtob(in); 
out = btopl(in); 

function [out] = mtob(in) 
% computing FMT from m to b.  
% in = m vector 
% out = b vector 
 
lm = length(in); 
natoms = round(log2(lm)); 		 
if 2^natoms == lm  
for step = 1:natoms  
	i124 = 2^(step-1); 			 
	i842 = 2^(natoms+1-step); 	 
	i421 = 2^(natoms - step); 	 
	in = reshape(in,i124,i842); 
	in(:,(1:i421)*2) = in(:,(1:i421)*2) + in(:,(1:i421)*2-1); 
end	 
out = reshape(in,1,lm); 
else 
	'ACCIDENT in mtob: length of input vector not OK: should be a power of 2' 
end 

function [out] = btopl(in) 
% compute pl from b 
% in = b 
% out = pl 
 
lm = length(in); 
natoms = round(log2(lm)); 		 
if 2^natoms == lm  
in = in(lm)-fliplr(in); 
in(1) = 0; 
out = in; 
else 
	'ACCIDENT in btopl: length of input vector not OK: should be a power of 2' 
end 


