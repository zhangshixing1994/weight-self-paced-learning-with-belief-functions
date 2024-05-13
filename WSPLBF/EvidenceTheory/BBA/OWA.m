
function [mass_f] = OWA(C,method)
%C是一个k*n的决策，每一列表示一个专家做出的一个决策

[k,n]=size(C);
cPess=min(C,[],2);
cOpti=max(C,[],2);
E=[cPess';cOpti'];

switch method
    case 1
        %cowa-er
        E_N=E/max(max(E));
        mass=zeros(k,2^k);
        for i=1:k
            mass(i,2^(i-1)+1)=E_N(1,i);
            mass(i,2^k-2^(i-1))=1-E_N(2,i);
            mass(i,2^k)=E_N(2,i)-E_N(1,i);
        end
        mass_f=DSFusionMatrix(mass,2);
    case 2
        %fcowa-er
        E(1,:)=E(1,:)/max(E(1,:));
        E(2,:)=E(2,:)/max(E(2,:));
        
        mpess=alpha_cut(E(1,:));
        mopti=alpha_cut(E(2,:));
        mass=[mpess;mopti];
        mass_f=DSFusionMatrix(mass,2);
end
    
end

function [m] = alpha_cut(v)
k=length(v);
v=v/max(v);
m=zeros(1,2^k);

w=[];
for i=1:k
    wi=2^(i-1);
    w=[w,wi];
end

w=w(1:k);
flag=zeros(1,k);

[value_sort,idx_sort]=sort(v,'descend');
value=value_sort(1,1:end-1)-value_sort(1,2:end);

for i=1:length(value)
    flag(idx_sort(i))=1;
    idx=sum(w.*flag)+1;
    m(1,idx)=value(i);
end
m(1,end)=value_sort(end);

end
        