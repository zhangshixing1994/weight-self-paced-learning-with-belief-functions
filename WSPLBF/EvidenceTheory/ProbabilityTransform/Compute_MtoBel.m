
function [Bel] = Compute_MtoBel(Mass)


%����������BetP
n=size(Mass,1);
Bel=[];
for i=1:n
bel=mtobel(Mass(i,:));
Bel=[Bel;bel];
end

end