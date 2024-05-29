function [ax] = plotRRM(matrix)

ax = imagesc(matrix);
line([.5 .5 6.5 6.5 .5],[.5 6.5 6.5 .5 .5], 'color', 'k', 'linewidth', 3);
line([6.5 6.5 12.5 12.5 6.5],[6.5 12.5 12.5 6.5 6.5], 'color', 'k', 'linewidth', 3);
%axis ij;
daspect([1 1 1]);
yticks(1:12);
xticks(1:12);
xlabel('trial b');
ylabel('trial a');
colormap jet;
% xticklabels(trList);
% yticklabels(trList);