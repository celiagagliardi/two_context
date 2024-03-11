%% Marc figures

panels = [{'Day 1'}, {'Day 2'}; {'Day 2'}, {'Day 3'}; {'Day 1'}, {'Day 3'}];
% replicate main figure
data = readtable('train_predict_different_days_cells_any_digs_Geo_Corr_kernel_polynomial_minCells_2.xlsx');
all_avg = nan(1,3); all_sem = nan(1,3);
figure;
for t = 1:size(panels,1)
    fprintf('All cells train %s and predict %s\n', panels{t,1}, panels{t,2});
    tbl = data(ismember(data.trainGroupLabel, panels(t,1)) & ismember(data.predictGroupLabel, panels(t,2)),:);
    overall_all = mean(tbl.accuracy);
    sem = std(tbl.accuracy) ./ sqrt(length(tbl.accuracy));
    all_avg(t) = overall_all;
    all_sem(t) = sem;
    [acc, sidx] = sort(tbl.accuracy, 'ascend');
    animals = tbl.animalName(sidx);
    subplot(1,3,t)
    barh([overall_all; acc]);
    hold on;
    xline(.5);
    yticklabels(['overall'; animals])
    errorbar(overall_all, 1, sem, 'horizontal', 'k');
    hold off;
    title(sprintf('Train %s, Predict %s', panels{t,1}, panels{t,2}));

    [h,p] = ttest(acc, .5, 'tail', 'right');
    fprintf('\t Overall Avg = %.2f, sem = %.2f \n', overall_all, sem);
    fprintf('\t\t ttest: h = %d, p = %.3f\n', h, p);
    [p,h] = signrank(acc, .5, 'tail', 'right');
    fprintf('\t\t sign rank: h = %d, p = %.3f\n', h, p);
    h = kstest(acc);
    if h
        fprintf('\t\t\tdata not normal, recommend sign rank test\n\n');
    end
end
sgtitle('Using all cells');
saveas(gcf, fullfile('figures', 'allcells_predictions.png'));

%%
data = readtable('train_predict_different_days_cells_stable_digs_Geo_Corr_kernel_polynomial_minCells_2.xlsx');
fi_avg = nan(1,3); fi_sem = nan(1,3);
figure;
for t = 1:size(panels,1)
    fprintf('FI cells train %s and predict %s\n', panels{t,1}, panels{t,2});
    tbl = data(ismember(data.trainGroupLabel, panels(t,1)) & ismember(data.predictGroupLabel, panels(t,2)),:);
    overall_FI = mean(tbl.accuracy);
    sem_FI = std(tbl.accuracy) ./ sqrt(length(tbl.accuracy));
    fi_avg(t) = overall_FI;
    fi_sem(t) = sem_FI;
    [acc, sidx] = sort(tbl.accuracy, 'ascend');
    animals = tbl.animalName(sidx);
    subplot(1,3,t)
    barh([overall_FI; acc]);
    hold on;
    xline(.5);
    yticklabels(['overall'; animals])
    errorbar(overall_FI, 1, sem_FI, 'horizontal', 'k');
    hold off;
    title(sprintf('Train %s, Predict %s', panels{t,1}, panels{t,2}));

    [h,p] = ttest(acc, .5, 'tail', 'right');
    fprintf('\t Overall Avg = %.2f, sem = %.2f \n', overall_FI, sem_FI);
    fprintf('\t\t ttest: h = %d, p = %.3f\n', h, p);
    [p,h] = signrank(acc, .5, 'tail', 'right');
    fprintf('\t\t sign rank: h = %d, p = %.3f\n', h, p);
    h = kstest(acc);
    if h
        fprintf('\t\t\tdata not normal, recommend sign rank test\n\n');
    end
end
sgtitle('Using FI cells');
saveas(gcf, fullfile('figures', 'ficells_predictions.png'));


%%
data = readtable('train_predict_different_days_cells_unstable_digs_Geo_Corr_kernel_polynomial_minCells_2.xlsx');
fs_avg = nan(1,3); fs_sem = nan(1,3);
figure;
for t = 1:size(panels,1)
    fprintf('FS cells train %s and predict %s\n', panels{t,1}, panels{t,2});
    tbl = data(ismember(data.trainGroupLabel, panels(t,1)) & ismember(data.predictGroupLabel, panels(t,2)),:);
    overall_FS = nanmean(tbl.accuracy);
    sem_FS = std(tbl.accuracy, 'omitnan') ./ sqrt(sum(~isnan(tbl.accuracy)));
    fs_avg(t) = overall_FS;
    fs_sem(t) = sem_FS;
    [acc, sidx] = sort(tbl.accuracy, 'ascend');
    animals = tbl.animalName(sidx);
    subplot(1,3,t)
    barh([overall_FS; acc]);
    hold on;
    xline(.5);
    yticklabels(['overall'; animals])
    errorbar(overall_FS, 1, sem_FS, 'horizontal', 'k');
    hold off;
    title(sprintf('Train %s, Predict %s', panels{t,1}, panels{t,2}));

    [h,p] = ttest(acc, .5, 'tail', 'right');
    fprintf('\t Overall Avg = %.2f, sem = %.2f \n', overall_FS, sem_FS);
    fprintf('\t\t ttest: h = %d, p = %.3f\n', h, p);
    [p,h] = signrank(acc, .5, 'tail', 'right');
    fprintf('\t\t sign rank: h = %d, p = %.3f\n', h, p);
    h = kstest(acc);
    if h
        fprintf('\t\t\tdata not normal, recommend sign rank test\n\n');
    else
        fprintf('\t\t\tdata normal, recommend ttest \n\n');
    end
end
sgtitle('Using FS cells');
saveas(gcf, fullfile('figures', 'fscells_predictions.png'));


%% Plot overall only
avg = [all_avg; fi_avg; fs_avg];
sem = [all_sem; fi_sem; fs_sem];



figure
hBar = bar(avg);
hBar(1).XData
hold on
for k1 = 1:size(avg,2)
    ctr(k1,:) = bsxfun(@plus, hBar(k1).XData, hBar(k1).XOffset');   % Note: ‘XOffset’ Is An Undocumented Feature, This Selects The ‘bar’ Centres
    ydt(k1,:) = hBar(k1).YData;                                     % Individual Bar Heights
end
hold on
errorbar(ctr, ydt, sem, 'k', 'linestyle', 'none')

set(gca, 'XTickLabel', {'All', 'FI', 'FS'});

title('Overall Heading prediction across time')
legend({'Train Day 1, Predict Day 2', ...
    'Train Day 2, Predict Day 3', ...
    'Train Day 1, Predict Day 3'});

saveas(gcf, fullfile('figures', 'overall_predictions.png'));


