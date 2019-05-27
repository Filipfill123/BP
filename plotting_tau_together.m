% marix_tau_plotted_together = [];
% arms_added_together = [];
% j = 1;
% for i = 0.610208:0.05:3
%     marix_tau_plotted_together(j) = critFcn(ksi*i, trajectory);
%     arms_added_together(j) = L1*i + L2*i + L3*i;
%     j = j + 1;
% end

tic
marix_tau_plotted_together = [];
arms_added_together = [];
j = 1;
for i = 0.610208:0.01:3
    marix_tau_plotted_together(j) = critFcn(ksi*i, trajectory);
    arms_added_together(j) = L1*i + L2*i + L3*i;
    j = j + 1;
end
 
%%
% Nalezení nejmenšího (a tedy nejlepšího) tau

smallest_tau_together = min(min(marix_tau_plotted_together));
[~, best_column] = find(marix_tau_plotted_together == smallest_tau_together);
best_L_1_together = arms_added_together(best_column)/3;
best_L_2_together = arms_added_together(best_column)/3;
best_L_3_together = arms_added_together(best_column)/3;
best_parametrs_and_tau_together = [best_L_1_together; best_L_2_together; best_L_3_together; smallest_tau_together]
toc

%%
% Vykreslení jenom do 2D - x = souèet ramen, y = tau
figure;
hold on;
scatter(arms_added_together, marix_tau_plotted_together,25);

scatter(arms_added_together(best_column), marix_tau_plotted_together(best_column),45,'r','filled');

title('Závislost tau na délce ramen');
xlabel('Souèet ramen L1 + L2 + L3');
ylabel('Tau');
grid;