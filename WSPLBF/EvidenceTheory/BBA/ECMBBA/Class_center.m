function [gplus,c] = Class_center(label_data,label_label)
%函数功能：找出每个焦元的重心

%输入：label_data,label_label已标记样本集与样本标签

%输出：gplus焦元重心(没有空集）
%输出：c每个焦元的元素的个数


K=length(unique(label_label));
[~,m]=size(label_data);
nbFoc=2^K;



%算出单点焦元重心g
ii=1:2^K;
F=zeros(length(ii),K);
for i=1:K,
      F(:,i)=bitget(ii'-1,i);
end;
c = sum(F(2:end,:),2)';  %焦元的个数
g =zeros(K,m);
for i=1:K
    index=find(label_label==i);
    g(i,:)=sum(label_data(index',:),1)/length(index);
end

%算出所有非空焦元重心gplus
gplus=[];
for i=2:nbFoc
    fi = F(i,:);
    truc = repmat(fi',1,m);
    gplus = [gplus;sum(g.*truc)./sum(truc)];%求得重心
end

end