function [gplus,c] = Class_center(label_data,label_label)
%�������ܣ��ҳ�ÿ����Ԫ������

%���룺label_data,label_label�ѱ����������������ǩ

%�����gplus��Ԫ����(û�пռ���
%�����cÿ����Ԫ��Ԫ�صĸ���


K=length(unique(label_label));
[~,m]=size(label_data);
nbFoc=2^K;



%������㽹Ԫ����g
ii=1:2^K;
F=zeros(length(ii),K);
for i=1:K,
      F(:,i)=bitget(ii'-1,i);
end;
c = sum(F(2:end,:),2)';  %��Ԫ�ĸ���
g =zeros(K,m);
for i=1:K
    index=find(label_label==i);
    g(i,:)=sum(label_data(index',:),1)/length(index);
end

%������зǿս�Ԫ����gplus
gplus=[];
for i=2:nbFoc
    fi = F(i,:);
    truc = repmat(fi',1,m);
    gplus = [gplus;sum(g.*truc)./sum(truc)];%�������
end

end