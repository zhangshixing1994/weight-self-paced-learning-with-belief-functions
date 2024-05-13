

function Ed= calc_ed(M,y)
%%�������ܣ�������������֮���֤�ݾ������

%���룺MΪ֤�ݷ������
%�����Evidence_distance_MatrixΪ��������֮���֤�ݾ���
[m,n]=size(M);
k=log2(n);
Dj= compute_Dj(k);
for i=1:m
    y_=zeros(1,n);
    y_(1,2^(y(i)-1)+1)=1;
    Ed(i)=sqrt((M(i,:)-y_)*Dj*(M(i,:)-y_)'*0.5);
 end
end



function D = compute_Dj(k)
%�������ܣ�����Jousselme������Ҫ�ľ���Dj
%���룺kΪ�������
%���������Jousselme������Ҫ�ľ���D

D=zeros(2^k,2^k);
for i=1:2^k
    for j=1:2^k
        Binary_Array1=zeros(1,k);%i�Ķ���������
        Binary_Array2=zeros(1,k);%j�Ķ���������
        for p=1:k
            Binary_Array1(1,p)=bitget(i-1,p);
            Binary_Array2(1,p)=bitget(j-1,p);
        end
        Binary_Array3=zeros(2,k);%��һ�б���i��j���룬�ڶ��б���i��j�Ļ�
        Binary_Array3(1,:)=Binary_Array1(1,:)&Binary_Array2(1,:);
        Binary_Array3(2,:)=Binary_Array1(1,:)|Binary_Array2(1,:);
        D(i,j)=sum(Binary_Array3(1,:))/sum(Binary_Array3(2,:));
    end
end
D(1,1)=0;%�ѵ�һ��Ԫ����0��
end