function [BetP] = M2BetP(mass)
[p,~]=size(mass);
BetP=[];
for i=1:p
[betp] = mtobetp2(mass(i,:));
BetP=[BetP;betp];
end

end

function [out] = mtobetp2(mm)
% computing BetP on W from the m vector 
% out = BetP vector: order a,b,c,... only on singletons
% I add a rescaling factor as with ccard >= 7 I had trouble beause m0 = 1.

%
%   Version : Nov 5,03.
%   Copyright (c) 2002 by Philippe Smets, Brussels, Belgium.
% 	Author: Philippe Smets
%
lm = length(mm);
natoms = round(log2(lm)); 

if mm(1) == 1
	out = 	ones(1,natoms)./natoms; % case of total conflict
else
	cf = [0 1];
	for i = 2:natoms
		cf = [cf cf+1];
	end
	for i = 1:natoms
		out(i) = (bitget(1:lm-1,i)./cf([2:lm]))*mm(2:lm)';
	end
	out = out/sum(out);
end
end

