% %%
% % sestrojení matice pro vykreslení
% matrix_tau_plotted_separately = [];
% pocet = 0;
% possible = 0;
% impossible = 0;
% column = 1;
% for i = 0.5:0.01:4
%     for j = 0.5:0.01:4
%         for z = 0.5:0.01:4
%             try
%                 matrix_tau_plotted_separately(1,column) = i;
%                 matrix_tau_plotted_separately(2,column) = j;
%                 matrix_tau_plotted_separately(3,column) = z;
%                 matrix_tau_plotted_separately(4,column) = critFcn([L1*i, L2*j, L3*z], trajectory);
%                 possible = possible + 1;
%                 probehlo = ['Úspìšnì probìhlo ' num2str(possible) ' výpoètù.'];
%                 disp(probehlo);
%                 
%             catch
%                 warning('Realsqrt produced complex result');
%                 matrix_tau_plotted_separately(4,column) = inf;
%                 impossible = impossible + 1;
%                 disp(impossible);
%             end    
%             column = column + 1;
%             pocet = pocet + 1;
%             disp(pocet);
%         end
%     end
% end
tic
matrix_tau_plotted_separately = [];
pocet = 0;
possible = 0;
impossible = 0;
column = 1;
for i = 0.02:0.05:1.2
    for j = 0.02:0.05:1.2
        for z = 0.02:0.05:1.2
            try
                matrix_tau_plotted_separately(1,column) = i;
                matrix_tau_plotted_separately(2,column) = j;
                matrix_tau_plotted_separately(3,column) = z;
                matrix_tau_plotted_separately(4,column) = critFcn([L1*i, L2*j, L3*z], trajectory);
                possible = possible + 1;
                probehlo = ['Úspìšnì probìhlo ' num2str(possible) ' výpoètù.'];
                disp(probehlo);
                
            catch
                warning('Realsqrt produced complex result');
                matrix_tau_plotted_separately(4,column) = inf;
                impossible = impossible + 1;
                disp(impossible);
            end    
            column = column + 1;
            pocet = pocet + 1;
            disp(pocet);
        end
    end
end

%%
% Nalezení nejmenšího tau

smallest_tau_separately = min(min(matrix_tau_plotted_separately(4,:)));
[~, best_column] = find(matrix_tau_plotted_separately(4,:) == smallest_tau_separately);
best_L_1_separately = matrix_tau_plotted_separately(1,best_column);
best_L_2_separately = matrix_tau_plotted_separately(2,best_column);
best_L_3_separately = matrix_tau_plotted_separately(3,best_column);
best_parametrs_and_tau_separately = [best_L_1_separately; best_L_2_separately; best_L_3_separately; smallest_tau_separately]
toc

%%
% Vykreslení
% 
% osaX = matrix_tau_plotted_separately(1,:);
% osaY = matrix_tau_plotted_separately(2,:);
% osaZ = matrix_tau_plotted_separately(3,:);
% osaTau = matrix_tau_plotted_separately(4,:);
% 
% scatter(osaX,osaY,osaZ,osaTau);
% grid;
% ax = gca;
% ax.XDir = 'reverse';
% view(-30,14)
% xlabel('Rameno L1');
% ylabel('Rameno L2');
% zlabel('Rameno L3');
% 
% cb = colorbar;
% cb.Label.String = 'Tau závislé na délce ramen';