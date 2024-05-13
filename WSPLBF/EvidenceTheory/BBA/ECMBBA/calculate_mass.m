function [Mass] = calculate_mass(x,gplus,c,Class_num,beta,alpha,delta)
%�������ܣ���������mass

%���룺x��������label������ǩ��Class_num����������
%���룺gplus ÿһ�������
%���룺c ��Ԫ����Ԫ�ظ������������ռ���

%�����Mass,BetP,AM



K=Class_num;
nbFoc=2^K;
[n,~]=size(x);
% 
% beta=2;alpha = 1;
% delta = 10 ;
delta2=delta^2;

%��������mass
D=pdist2(x,gplus).^(-2/(beta-1));
C=repmat(c.^((-alpha/(beta-1))),n,1);
CC=repmat(sum(D.*C,2)+delta2^(-beta+1),1,nbFoc-1);
m=D.*C./ CC;
m_empty=1-sum(m,2);
Mass=[m_empty,m];




end



