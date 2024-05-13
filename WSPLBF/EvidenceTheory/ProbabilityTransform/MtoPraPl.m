


function [PraPl] = MtoPraPl(mass,Class_num)

Pl=Compute_MtoPl(mass);
Bel=Compute_MtoBel(mass);


[m,~]=size(mass);
flag=fliplr(binaryMatrix(Class_num));
PraPl=[];

singlton_idx=[];
for k=1:Class_num
    singlton_idx=[singlton_idx,2^(k-1)+1];
end

PraPl=[];
for i=1:m

    epsilon=(1-sum(Bel(i,singlton_idx)))/(sum(Pl(i,singlton_idx)));
    prapl_i=Bel(i,singlton_idx)+ epsilon*Pl(i,singlton_idx);
    PraPl=[PraPl;prapl_i];
end
v=sum(PraPl,2);
D=diag(v);
PraPl=D^-1*PraPl;

end

