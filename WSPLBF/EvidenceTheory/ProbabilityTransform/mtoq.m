function [out] = mtoq(in)
% computing FMT from m to q. 
% q是一个中间变量，1*lm的矩阵，第一个元素存储所有包含空集的焦元的mass函数之和
% q(A)=m(B1)+m(B2)+...+m(Bn),A是Bi的子集
% in = m vector
% out = q vector

lm = length(in);
natoms = round(log2(lm)); 		
if 2^natoms == lm 
for step = 1:natoms 
	i124 = 2^(step-1); 			
	i842 = 2^(natoms+1-step); 	
	i421 = 2^(natoms - step); 	
	in = reshape(in,i124,i842);
	in(:,(1:i421)*2-1) = in(:,(1:i421)*2-1) + in(:,(1:i421)*2);
end	
out = reshape(in,1,lm);
else
	'ACCIDENT in mtoq: length of input vector not OK: should be a power of 2'
end
