function [w] = calc_w(betp,y,beta)


[n,k]=size(betp);
w=[];
for i=1:n
    w=[w;betp(i,y(i))];
end
w=w.^beta;
end

