function [AM,w] = calc_evidinfo(Xtrain,ytrain)

k=5;
mass= EKNNBBA_train(Xtrain,ytrain,k);
[AM,BetP] = compute_AM(mass,ytrain);
gama=1;
[w] = calc_w(BetP,ytrain,gama);

end

