function [idx_train,idx_test] =train_test_split2(datalabel,sample_rate)

k=length(unique(datalabel));%类别个数

idx_train=[];
idx_test=[];
for i=1:k

    index=find(datalabel==i);%找到第i个类的索引
    index=index(randperm(numel(index)));%打乱索引
    train_idx=index(1:ceil(length(index)*sample_rate));%取索引的前半部分为训练索引
    test_idx=index(ceil(length(index)*sample_rate)+1:length(index));%取索引的后半部分为测试索引
    idx_train=[idx_train;train_idx];
    idx_test=[idx_test;test_idx];
end

end

