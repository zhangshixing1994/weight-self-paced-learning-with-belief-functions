
function [bestParam] = search_para(Xtrain,ytrain,Xtest,ytest,net,AM,xi,w)



% 
gama= optimizableVariable('gama',  [0 100], 'Type', 'real');
psi=optimizableVariable('psi',  [0 10], 'Type', 'real');


parameter = [gama,psi];

% 目标函数
objFun = @(parameter) fun(parameter,Xtrain,ytrain,Xtest,ytest,net,AM,xi,w);

% 贝叶斯优化
iter = 30;
points = 10;
results = bayesopt(objFun, parameter, 'Verbose', 0, ...
                   'MaxObjectiveEvaluations', iter,...
                   'NumSeedPoints', points,'PlotFcn',[]);%
               % 优化结果
[bestParam, ~, ~] = bestPoint(results, 'Criterion', 'min-observed');
end

function [objective] = fun(parameter,Xtrain,ytrain,Xtest,ytest,net,AM,xi,w)


gama= parameter.gama;
psi=parameter.psi;
w=w.^gama;
xi=xi.^psi;

output = WSPLBF(Xtrain,ytrain,Xtest,ytest,net,AM,xi,w);

acc=output.acc;

objective =1-acc;

end
