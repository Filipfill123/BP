fun = @critFcn_2;
x0 = [0.1 1.5 3];
lb = [0.01 0.01 0.01];
ub = [2 2 2];
options = optimoptions('patternsearch','PlotFcn',@psplotbestf);

pattern_search_tStart = tic;
[ksi_pattern_search, fval_pattern_search,exitFlag,Output_pattern_search] = patternsearch(fun,x0,[],[],[],[],lb,ub,[],options);
pattern_search_tEnd = toc(pattern_search_tStart);
%%
fprintf('The number of iterations was : %d\n', Output_pattern_search.iterations);
pattern_search_funccount = Output_pattern_search.funccount;
fprintf('The number of function evaluations was : %d\n', pattern_search_funccount);
fprintf('The best function value found was : %g\n', fval_pattern_search);
g=sprintf('%d ', ksi_pattern_search);
fprintf('The best ksi found by simulated annealing was: %s\n', g);