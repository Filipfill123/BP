time_ellapsed = [gen_alg_tEnd,part_swarm_tEnd,sim_ann_tEnd,pattern_search_tEnd];
c = categorical({'GA','PSO','SA','PSA'});
bar(c,time_ellapsed,0.4);
ylim([100 900])
title('Èas potøebný k probìhnutí algoritmu');