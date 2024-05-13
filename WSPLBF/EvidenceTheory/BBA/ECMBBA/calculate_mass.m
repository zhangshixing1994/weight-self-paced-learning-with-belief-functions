function [Mass] = calculate_mass(x,gplus,c,Class_num,beta,alpha,delta)
%函数功能：计算样本mass

%输入：x样本集，label样本标签，Class_num样本类别个数
%输入：gplus 每一类的重心
%输入：c 焦元包含元素个数（不包含空集）

%输出：Mass,BetP,AM



K=Class_num;
nbFoc=2^K;
[n,~]=size(x);
% 
% beta=2;alpha = 1;
% delta = 10 ;
delta2=delta^2;

%计算样本mass
D=pdist2(x,gplus).^(-2/(beta-1));
C=repmat(c.^((-alpha/(beta-1))),n,1);
CC=repmat(sum(D.*C,2)+delta2^(-beta+1),1,nbFoc-1);
m=D.*C./ CC;
m_empty=1-sum(m,2);
Mass=[m_empty,m];




end



