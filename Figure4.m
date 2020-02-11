%
% Filename    : Figure4.m
% Author      : Joseph Gleason
% Description : Generate Figure 4 for Automatica paper "Lagrangian
%               Approximations for Stochastic Reachability of a Target
%               Tube"               
%                       
%               Scalability Figure
% 

clearvars;
close all;

cmptimes = load('var/mats/HSCC_scalability_comptimes_02132019.mat');

MAXLENGTH = length(cmptimes.lagover.comptimes);

markers = struct('lagunder', '-g^', 'lagover', '-rd', ...
                 'dp', '-b>', 'genzps', '-cs', 'ccc', '-mo');
display_names = struct('lagunder', 'Lagrangian Underapproximation', ...
                       'lagover', 'Lagrangian Overapproximation', ...
                       'dp', 'Dynamic Programming', ...
                       'genzps', 'Genz + Patternsearch', ...
                       'ccc', 'Convex Chance Constrained');
             
flds = fields(cmptimes);
for lv = 1:length(flds)
    field_name = flds{lv};
    y = cmptimes.(field_name).comptimes;
    t = 2:min(length(y), MAXLENGTH)+1;
    y = y(1:length(t));
    plt = semilogy(t, y, markers.(field_name), ...
        'MarkerFaceColor', markers.(field_name)(2), ...
        'DisplayName', display_names.(field_name));
    
    if lv == 1
        hold on;
    end
end

hl = legend;
hl.Position = [0.1350 0.0135 0.7076 0.2007];
hl.FontSize = 8;
hl.Interpreter = 'latex';
hl.Box = 'off';

ha = gca;
ylim = ha.YLim;

yyaxis right;

ha.GridColor = 'k';
ha.MinorGridColor = 'k';
ha.XColor = 'k';
ha.YColor = 'k';

ha.YLim = ylim;
ha.FontSize = 9;

plot([2, MAXLENGTH+1], [1, 1], '--k', 'HandleVisibility', 'off')
plot([2, MAXLENGTH+1], [60, 60], '--k', 'HandleVisibility', 'off')
plot([2, MAXLENGTH+1], [3600, 3600], '--k', 'HandleVisibility', 'off')
plot([2, MAXLENGTH+1], [86400, 86400], '--k', 'HandleVisibility', 'off')

ha.YScale = 'log';
yticks([1, 60, 3600, 86400])
yticklabels({'', '1 min.', '1 hr.', '1 day'});
hold off;

ha.TickLabelInterpreter = 'latex';

yyaxis left;
ha.YLabel.Interpreter = 'latex';
ha.XLabel.Interpreter = 'latex';
ha.YLabel.Interpreter = 'latex';
ylabel('Computation Time [s]');
xlabel('Dimension')

hf = gcf;
hf.Units = 'inches';
hf.Position = [5.0729 2.7083 3 3.8438];

ha.Position = [0.1350 0.3268 0.7357 0.5982];

print(hf, '-r300', '-dpng', sprintf('./var/figs/figure4_%s.png', ...
    datestr(now, 'dd-mm-yyyy-HHMMSS')));