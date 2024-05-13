
function [Pl] =Compute_MtoPl(Mass)


%计算样本的BetP
n=size(Mass,1);
Pl=[];
for i=1:n
pl=mtopl(Mass(i,:));
Pl=[Pl;pl];
end

end