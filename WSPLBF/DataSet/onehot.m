function [y_] = onehot(y)
n=length(y);
k=length(unique(y));
y_=zeros(n,k);
for i=1:n
    y_(i,y(i))=1;
end

end

