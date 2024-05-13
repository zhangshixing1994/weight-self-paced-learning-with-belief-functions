

function [Diss] = compute_Diss(mass)


[pl] = Compute_MtoPl(mass);
%pl(find(pl==0))=10^-16;


%¼ÆËãÑù±¾Diss
Diss=-sum(mass(:,2:end).*log(pl(:,2:end))/log(2),2);



end