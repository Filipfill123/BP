%% Particle swarm optimization
%left boundary [0.01 0.01 0.01], right boundary [2 2 2], 3 variables
rng default % for reproducibility
ObjectiveFunction = @critFcn_2;
options = optimoptions(@particleswarm,'PlotFcn',@pswplotbestf);
%options = optimoptions(options, 'Display','iter');
lb = [0.01 0.01 0.01];
ub = [2 2 2];
nvars = 3;
part_swarm_tStart = tic;
[ksi_part_swarm,fval_part_swarm,exitFlag_part_swarm,Output_part_swarm] = particleswarm(ObjectiveFunction,nvars,lb,ub,options);
part_swarm_tEnd = toc(part_swarm_tStart);
%%
formatstring = 'Particleswarm reached the value %f using %d function evaluations.\n';
part_swarm_funccount = Output_part_swarm.funccount;
fprintf(formatstring,fval_part_swarm,part_swarm_funccount)
g=sprintf('%d ', ksi_part_swarm);
fprintf('The best ksi found by simulated annealing was: %s\n', g)
