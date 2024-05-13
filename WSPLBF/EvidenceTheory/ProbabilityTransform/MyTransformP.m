

function [MTP] = MyTransformP(Mass,prop)


%我的概率转换
n=size(Mass,1);
MTP=[];

Class_num=log2(length(Mass(1,:)));

flag=fliplr(binaryMatrix(Class_num));
for i=1:n
    mtp=[];
    for k=1:Class_num
        
        idx=find(flag(:,k)==1);
        mtp_k=sum(Mass(i,idx)*prop(idx,k));
        mtp=[mtp,mtp_k];
    end
MTP=[MTP;mtp];
    
    
end

end

