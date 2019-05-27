%% Simulated annealing ("simulované žíhání")
% defaultnì je teplota nastaveno na 100 pro každou dimenzi
% startingPoint nastaven na [0.1 0.1 0.1], left boundary na [0.01 0.01
% 0.01] a right boundary na [2 2 2]
rng default % For reproducibility
ObjectiveFunction = @critFcn_2;
startingPoint = [0.1 0.1 0.1];
lb = [0.01 0.01 0.01];
ub = [2 2 2];
options = optimoptions(@simulannealbnd, ...
                     'PlotFcn',{@saplotbestf,@saplottemperature,@saplotf,@saplotstopping});
options.ReannealInterval = 100;                 
sim_ann_tStart = tic;                 
[ksi_sim_ann,fval_sim_ann,exitFlag_sim_ann,Output_sim_ann] = simulannealbnd(ObjectiveFunction,startingPoint,lb,ub,options);
sim_ann_tEnd = toc(sim_ann_tStart);


fprintf('The number of iterations was : %d\n', Output_sim_ann.iterations);
sim_ann_funccount = Output_sim_ann.funccount;
fprintf('The number of function evaluations was : %d\n', sim_ann_funccount);
fprintf('The best function value found was : %g\n', fval_sim_ann);
g=sprintf('%d ', ksi_sim_ann);
fprintf('The best ksi found by simulated annealing was: %s\n', g);
