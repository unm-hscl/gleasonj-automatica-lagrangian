%
%

clearvars;
close all;

% load runtimes
rt = load('../examples/scalability-runtimes.mat');

% plot
semilogy(2:length(rt.under_runtimes)+1, rt.under_runtimes, '-kx', ...
    'MarkerFaceColor', 'k', 'MarkerSize', 7);
hold on;
semilogy(2:length(rt.under_runtimes)+1, rt.over_runtimes, '-ko', ...
    'MarkerFaceColor', 'k', 'MarkerSize', 4);
hold off;

xlabel('Dimension, $n$')
ylabel('Computation Time [s]')


% Formatting
FONT_SIZE = 8;


hf = gcf;
hf.Units = 'inches';
hf.Position(3) = 3;
hf.Position(4) = 1.5;

ha = gca;
ha.FontSize = FONT_SIZE;
ha.TickLabelInterpreter = 'latex';

ha.XLabel.FontSize = FONT_SIZE;
ha.XLabel.Interpreter = 'latex';

ha.YLabel.FontSize = FONT_SIZE;
ha.YLabel.Interpreter = 'latex';

ha.YTick = [1, 10, 100, 1000, 10000, 100000];
ha.XTick = [2, 5, 8, 11, 14];
ha.XLim = [2, 14];

hleg = legend({'$\hbox{Reach}_{0}^{\flat-}$', '$\hbox{Reach}_{0}^{\sharp+}$'}, ...
    'Interpreter', 'latex');
hleg.Position = [0.5801    0.2636    0.2845    0.1713];

print('-dpng', '-r450', 'coi-scalability.png');