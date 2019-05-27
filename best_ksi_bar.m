best_ksi = [ksi_gen_alg;ksi_part_swarm;ksi_sim_ann;ksi_pattern_search];
c = categorical({'GA','PSO','SA','PSA'});
bar(c,best_ksi,0.4);
ylim([0.001 1.2])
title('Nejlepší nalezená hodnota optimalizovaných parametrù');