

function [PrBl] = MtoPrBl(mass,Class_num)

Bel=Compute_MtoBel(mass);
[m,~]=size(mass);
flag=fliplr(binaryMatrix(Class_num));
PrBl=[];
for i=1:m
    prbl=[];
    for j=1:Class_num
        idx=find(flag(:,j)==1);
        prop=repmat(Bel(i,2^(j-1)+1),1,length(Bel(i,idx)))./Bel(i,idx);%±ÈÀý
        prbl_i=sum(prop.*mass(i,idx));
        prbl=[prbl,prbl_i];
    
    end
    PrBl=[PrBl;prbl];
end
v=sum(PrBl,2);
D=diag(v);
PrPl=D^-1*PrBl;

end
