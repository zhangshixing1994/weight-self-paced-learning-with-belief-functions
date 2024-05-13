

function Ed= calc_ed(M,y)
%%函数功能：计算样本两两之间的证据距离矩阵

%输入：M为证据分配矩阵
%输出：Evidence_distance_Matrix为样本两两之间的证据距离
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
%函数功能：计算Jousselme距离需要的矩阵Dj
%输入：k为聚类簇数
%输出：计算Jousselme距离需要的矩阵D

D=zeros(2^k,2^k);
for i=1:2^k
    for j=1:2^k
        Binary_Array1=zeros(1,k);%i的二进制数组
        Binary_Array2=zeros(1,k);%j的二进制数组
        for p=1:k
            Binary_Array1(1,p)=bitget(i-1,p);
            Binary_Array2(1,p)=bitget(j-1,p);
        end
        Binary_Array3=zeros(2,k);%第一行保存i和j的与，第二行保存i和j的或
        Binary_Array3(1,:)=Binary_Array1(1,:)&Binary_Array2(1,:);
        Binary_Array3(2,:)=Binary_Array1(1,:)|Binary_Array2(1,:);
        D(i,j)=sum(Binary_Array3(1,:))/sum(Binary_Array3(2,:));
    end
end
D(1,1)=0;%把第一个元素置0；
end