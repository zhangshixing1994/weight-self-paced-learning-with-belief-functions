function [data] = data_standard( data)
%���ݹ�һ��
[m,n]=size(data);
data=data-repmat(min(data),m,1);
data=data./repmat(max(data)-min(data),m,1);

end

