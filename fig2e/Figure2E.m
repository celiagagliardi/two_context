%% Figure 2E

data = readtable('train_predict_Geo_Corr_same_day_using_centeroutangledeg_nostability_accuracies.xlsx');
fprintf('Predicting Corr and Geo digs\n');
for d = 1:3
    dtbl = data(data.groupId ==d,:);
    cdata = dtbl.accuracy;
    avg = mean(cdata);
    sem = std(cdata) ./ sqrt(length(cdata));
    [h,p,ci, stats] = ttest(cdata, .5, 'tail', 'right');
    fprintf('Day %d: avg = %.3f, sem = %.3f, t(%d) = %.2f, p = %.3f\n', d, avg, sem, stats.df, stats.tstat,p);
    if kstest(cdata)
         fprintf('Data not normal, report sign rank test\n')
     p = signrank(cdata, .5, 'tail', 'right');
    fprintf('sign rank p = %.3f\n', p);
    end
end

%%

data = readtable('train_predict_Wrong_Feat_same_day_using_centeroutangledeg_nostability_accuracies.xlsx')
fprintf('Predicting Wrong and Feat digs\n');
for d = 1:3
    dtbl = data(data.groupId ==d,:);
    cdata = dtbl.accuracy;
    avg = mean(cdata);
    sem = std(cdata) ./ sqrt(length(cdata));
    [h,p,ci, stats] = ttest(cdata, .5, 'tail', 'right');
    fprintf('Day %d:avg = %.3f, sem = %.3f, one sample t test t(%d) = %.2f, p = %.3f\n', d, avg, sem, stats.df, stats.tstat,p);
    if kstest(cdata)
        fprintf('Data not normal, report sign rank test\n')
        p = signrank(cdata, .5, 'tail', 'right');
        fprintf('sign rank p = %.3f\n', p);
    end
end
