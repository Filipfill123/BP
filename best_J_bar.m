best_J = [fval_gen_alg,fval_part_swarm,fval_sim_ann,fval_pattern_search];
c = categorical({'GA','PSO','SA','PSA'});
bar(c,best_J,0.4);
ylim([99.6 99.8])
title('Nejlepší nalezená hodnota kriteriální funkce');

