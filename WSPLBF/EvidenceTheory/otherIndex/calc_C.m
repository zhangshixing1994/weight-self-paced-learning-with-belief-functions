function [C] = calc_C(k)
%�����ͻ����Ҫ�ľ���
%Ckl=1   if  Fk �� Fl�Ľ���Ϊ��
%Ckl=0       ���������
C=zeros(2^k,2^k);
for i=0:2^k-1
    for j=0:2^k-1
        result=bitand(i,j);%��λ����
        if(result==0)%�������Ϊ������1��������0
            C(i+1,j+1)=1;
        end
    end
end
end

