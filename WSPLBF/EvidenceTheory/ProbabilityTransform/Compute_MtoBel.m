
function [Bel] = Compute_MtoBel(Mass)


%计算样本的BetP
n=size(Mass,1);
Bel=[];
for i=1:n
bel=mtobel(Mass(i,:));
Bel=[Bel;bel];
end

end