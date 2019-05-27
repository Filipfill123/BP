function_eval = [gen_alg_funccount, part_swarm_funccount, sim_ann_funccount,pattern_search_funccount];
c = categorical({'GA','PSO','SA','PSA'});
bar(c,function_eval,0.4)
ylim([500 2650]);
title('Poèet vyhodnocení kriteriální funkce');