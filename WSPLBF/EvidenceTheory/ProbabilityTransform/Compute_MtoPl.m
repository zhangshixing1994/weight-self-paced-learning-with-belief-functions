
function [Pl] =Compute_MtoPl(Mass)


%����������BetP
n=size(Mass,1);
Pl=[];
for i=1:n
pl=mtopl(Mass(i,:));
Pl=[Pl;pl];
end

end