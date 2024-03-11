clear all; clc;
panels = {'Day 1', 'Day 2'; 'Day 1', 'Day 3'; 'Day 2', 'Day 3'};
data = readtable("calcium_stability_classification_FINAL-for Figure 6B.xlsx");
for t = 1:size(panels,1)
    tbl = data(ismember(data.dayA, panels{t,1}) & ismember(data.dayB, panels{t,2}),:);
    X = [tbl.x_stable_stable, tbl.x_stable_unstable];
    figure;
    subplot(1,2,1)
    pie(X);
    title('FI');
    legend('Remained FI', 'Became FS');

    X = [tbl.x_unstable_Unstable, tbl.x_unstable_stable];
    subplot(1,2,2)
    pie(X);
    title('FS');
    legend('Remained FS', 'Became FI');

    sgtitle(sprintf('%s to %s', panels{t,1}, panels{t,2}));
end




%% Figure 6C
data = readtable('train_predict_different_days_cells_stable_digs_Geo_Corr_kernel_polynomial_minCells_2.xlsx');


avg = nan(1,3);
sem = nan(1,3);
cData = cell(1,3);
for t = 1:size(panels,1)
    tbl = data(ismember(data.trainGroupLabel, panels{t,1}) & ismember(data.predictGroupLabel, panels{t,2}),:);
    acc = tbl.accuracy;
    cData{t} = acc;
    davg = mean(acc);
    dsem = std(acc) ./ sqrt(length(acc));
    avg(t) = davg;
    sem(t) = dsem;
    [h,p,ci, stats] = ttest(acc, .5, 'tail', 'right');
    fprintf('Train %s Predict %s\n', panels{t,1}, panels{t,2});
    fprintf('\t avg = %.3f, sem = %.3f\n', davg, dsem);
    fprintf('\t\t t test: h = %d, p = %.3f,  t(%d) = %.3f\n', h, p, stats.df, stats.tstat);
    [p,h] = signrank(acc, .5, 'tail', 'right');
    fprintf('\t\t sign rank: h = %d, p = %.3f\n', h, p);
    h = kstest(acc);
    if h
        fprintf('\t\t\t data not normal, recommend sign rank test\n');
    else
        fprintf('\t\t\t data normal, recommend t test\n');
    end
end

figure;
bar(avg);
hold on;
errorbar(avg, sem, 'k', 'linestyle', 'none');
for t = 1:3
    scatter(t, cData{t}, 'o', 'jitter', 'on');
end
yline(.5, 'k--');
xticklabels({sprintf('Train %s Predict %s', panels{1,1}, panels{1,2}), ...
    sprintf('Train %s Predict %s', panels{2,1}, panels{2,2}), ...
    sprintf('Train %s Predict %s', panels{3,1}, panels{3,2})});
ylabel('Prediction Accuracy');
