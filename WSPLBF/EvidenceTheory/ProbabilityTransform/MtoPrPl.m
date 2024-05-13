

function [PrPl] = MtoPrPl(mass,Class_num)

Pl=Compute_MtoPl(mass);
[m,~]=size(mass);
flag=fliplr(binaryMatrix(Class_num));
PrPl=[];
for i=1:m
    prpl=[];
    for j=1:Class_num
        idx=find(flag(:,j)==1);
        prop=repmat(Pl(i,2^(j-1)+1),1,length(Pl(i,idx)))./Pl(i,idx);%±ÈÀý
        prpl_i=sum(prop.*mass(i,idx));
        prpl=[prpl,prpl_i];
    
    end
    PrPl=[PrPl;prpl];
end
v=sum(PrPl,2);
D=diag(v);
PrPl=D^-1*PrPl;

end

