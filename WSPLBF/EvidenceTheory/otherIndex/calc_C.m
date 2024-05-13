function [C] = calc_C(k)
%计算冲突所需要的矩阵
%Ckl=1   if  Fk 与 Fl的交集为空
%Ckl=0       其他的情况
C=zeros(2^k,2^k);
for i=0:2^k-1
    for j=0:2^k-1
        result=bitand(i,j);%按位相与
        if(result==0)%如果交集为空则置1，否则置0
            C(i+1,j+1)=1;
        end
    end
end
end

