smallest_tau_separately = min(min(maticeTauPlottedSeparately(4,:)));
[~, best_column] = find(maticeTauPlottedSeparately(4,:) == smallest_tau_separately);
best_L_1_separately = maticeTauPlottedSeparately(1,best_column);
best_L_2_separately = maticeTauPlottedSeparately(2,best_column);
best_L_3_separately = maticeTauPlottedSeparately(3,best_column);
best_parametrs_and_tau_separately = [best_L_1_separately; best_L_2_separately; best_L_3_separately; smallest_tau_separately]