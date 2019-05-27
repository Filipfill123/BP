smallest_tau_together = min(min(maticeTauPlottedTogether));
[~, sloupec] = find(maticeTauPlottedTogether == smallest_tau_together);
best_L_1_together = soucetRamen(sloupec)/3;
best_L_2_together = soucetRamen(sloupec)/3;
best_L_3_together = soucetRamen(sloupec)/3;
best_parametrs_and_tau_together = [best_L_1_together; best_L_2_together; best_L_3_together; smallest_tau_together]