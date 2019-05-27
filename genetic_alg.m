
%% Genetický algoritmus 
% population nastavena automaticky na 50, MaxGeneretions empiricky
% nastaveno na 38, MaxStallGenerations nastaveno na 15
%3 variables, left boundary [0.05 0.05 0.05] a right boundary [4 4 4]
rng default % For reproducibility
FitnessFunction = @critFcn_2;
%'Display','iter'
opts = optimoptions(@ga,'PlotFcn',{@gaplotbestf,@gaplotstopping},'Display','iter');
opts = optimoptions(opts,'MaxGenerations',50,'MaxStallGenerations', 15);

numberOfVariables = 3;
lb = [0.01,0.01,0.01];
ub = [2,2,2];
gen_alg_tStart = tic;
[ksi_gen_alg,fval_gen_alg,exitFlag_gen_alg,Output_gen_alg] = ga(FitnessFunction,numberOfVariables,[],[],[],[],lb,ub,[],opts);
gen_alg_tEnd = toc(gen_alg_tStart);

fprintf('The number of generations was : %d\n', Output_gen_alg.generations);
fprintf('The best function value found was : %g\n', fval_gen_alg);
g=sprintf('%d ', ksi_gen_alg);
fprintf('The best ksi found by genetic algorithm was: %s\n', g);
%%
gen_alg_funccount = Output_gen_alg.funccount;
fprintf('The number of function evaluations was : %d\n', Output_gen_alg.funccount);




  
  