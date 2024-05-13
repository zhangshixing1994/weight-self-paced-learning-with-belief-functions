function [idx_train,idx_test] =train_test_split2(datalabel,sample_rate)

k=length(unique(datalabel));%������

idx_train=[];
idx_test=[];
for i=1:k

    index=find(datalabel==i);%�ҵ���i���������
    index=index(randperm(numel(index)));%��������
    train_idx=index(1:ceil(length(index)*sample_rate));%ȡ������ǰ�벿��Ϊѵ������
    test_idx=index(ceil(length(index)*sample_rate)+1:length(index));%ȡ�����ĺ�벿��Ϊ��������
    idx_train=[idx_train;train_idx];
    idx_test=[idx_test;test_idx];
end

end

