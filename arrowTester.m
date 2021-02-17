close all;
figure('Position', [0, 0, 1440, 760])
set(gca, 'Position', [0 0 1 1])
xlim([-5, 5]);
ylim([-5, 5]);

% % myArrow([0, 1], [0, 1], 1.5, 50, 50, 'k', '-', 0.1*pi)
% for i = -pi:pi/10:pi
%     cla
%     myArrow([0, 1], [0, 0], 1.5, 10, 10, 'k', '-', i, false)
%     drawnow;
% end
    myArrow([0, 1], [0, 0], 1.5, 10, 10, 'k', '-', -0.5*pi, false)

% annotation('arrow', [0.5, 0.5], [0.6, 0.6])