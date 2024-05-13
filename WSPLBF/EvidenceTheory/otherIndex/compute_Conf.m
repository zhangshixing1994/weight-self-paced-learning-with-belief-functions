
function [Conf] = compute_Conf(mass)


[bel] = Compute_MtoBel(mass);
bel(find(bel==0))=10^-16;


%¼ÆËãÑù±¾Conf
Conf=-sum(mass(:,2:end).*log(bel(:,2:end))/log(2),2);

end