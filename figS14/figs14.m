%% FigS14
% This script generates the supplemental figure that complements figure 6,
% showing the change in cell identity as a function of correlation
% threshold (that determines the cell identity).

% FS cells
data = readtable('figure_6B_alternate_cellreg_count_results_FS_remained_FS.xlsx');

figure;
for p = 1:size(data,1)
    subplot(10,3,p)
    pie([data.pChanged(p), data.pRemained(p)])
    title({sprintf('Decile %d', data.BIN_ID_A(p)), sprintf('%s to %s', data.groupLabelA{p}, data.groupLabelB{p})})
end

% FI cells
data = readtable('figure_6B_alternate_cellreg_count_results_FI_remained_FI.xlsx');
figure;
for p = 1:size(data,1)
    subplot(10,3,p)
    pie([data.pChanged(p), data.pRemained(p)])
    title({sprintf('Decile %d', data.BIN_ID_A(p)), sprintf('%s to %s', data.groupLabelA{p}, data.groupLabelB{p})})
end

%%
reg = load('calcium_cellregscores_2plusdays_per_decile_20240129.csv');
figure;
boxplot(reg)


%% get color scheme
a = [1 2 3 4 5 6 7 8 9 10];
figure
imagesc(a);
colormap jet;

pcolor(a)