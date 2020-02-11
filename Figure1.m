%
% Name: Figure1.m
% Author: Joseph Gleason
% Description: Draw the set expansion figure.
%

% meshgrid
[X, Y] = meshgrid(-1.3:0.005:1.3, -1.3:0.005:1.3);
Z1 = mvnpdf([X(:), Y(:)], [0.3, 0], 0.06*eye(2));

Z2 = mvnpdf([X(:), Y(:)], [-0.5, 0], 0.04*eye(2));
P1 = Polyhedron([-0.8, 0; -0.3, -0.6; 0.3, 0; 0, 0.4; -0.4, 0.4]);
P2 = 1.4 * P1;
P1z = Polyhedron(horzcat(P1.V, [0.001; 0.001; 0.001; 0.001; 0.001]));

C = zeros([size(X), 3]);
for i = 1:numel(X)
    v = [X(i); Y(i)];
    row = mod(i, size(X, 1));
    if row == 0; row = size(X, 1); end
    col = floor((i-1) / size(X, 1)) + 1;
    if P1.contains(v)
        C(row, col, :) = [1, 0, 0];
%     elseif P2.contains(v)
%         C(row, col, :) = [0, 0, 1];
    else
        C(row, col, :) = [1, 1, 0];
    end
end

Z = Z1 + Z2;
Z = reshape(Z, size(X));

sf = surf(X, Y, Z, C);
sf.EdgeColor = 'none'; 
% sf.FaceColor = [0, 0.7, 0];
sf.FaceAlpha = 0.5;
sf.FaceLighting = 'gouraud';
sf.AlphaData = 0.5 * ones(size(Z));
% sf.FaceAlpha = 'interp';

light(gca, 'Position', [-1, -1, 10], 'Style', 'infinite', 'Color', [1, 1, 1]);

hold on;
% patch(1.5 * [-1, 1, 1, -1], 1.5* [-1, -1, 1, 1], 0.1 * [-1, -1, -1, -1], [1, 1, 0]);
plot(P2, 'Color', [0, 0, 1])
plot(P1z)
hold off;

grid off;
ha = gca;
ha.TickLength = [0, 0];
ha.XTickLabel = '';
ha.YTickLabel = '';
ha.ZTickLabel = '';

%% Formatting for paper
hf = gcf;
hf.Units = 'inches';
hf.Position = [hf.Position(1:2), 3, 3];

ha.FontSize = 9;
ha.XLabel.Interpreter = 'latex';
ha.XLabel.String = '$x_{1}$';
ha.YLabel.Interpreter = 'latex';
ha.YLabel.String = '$x_{2}$';
ha.ZAxis.Visible = 'off';

Pt = annotation(hf, 'textbox', 'String', '$\mathcal{P}$', ...
    'Interpreter', 'latex', 'FontSize', 9, 'LineStyle', 'none', ...
    'Position', [0.4556 0.2403 0.1019 0.0847]);

TPt = annotation(hf, 'textbox', 'String', '$\theta^{*}\mathcal{P}$', ...
    'Interpreter', 'latex', 'FontSize', 9, 'LineStyle', 'none', ...
    'Position', [0.4056 0.1674 0.1019 0.0847]);

