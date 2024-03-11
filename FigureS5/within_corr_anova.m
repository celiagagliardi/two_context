%% within_corr anova
data = readtable('caltet_pairwise_correlation_scatters_by_bestaligned_data_20230615_104325.xlsx');
decs = discretize(data.bestAlignedCorrelation, quantile(data.bestAlignedCorrelation, [0:10]/10));

data(ismember(decs,1),:);

T = [data, table(decs, 'variableName', {'decile'})];
writetable(T, 'across_within_corr_figureS5.xlsx');

[p,anovatab,stats] = anova1(T.avgWithinCorrelation, T.decile);
anovaResults = multcompare(stats)

results = array2table(anovaResults, 'variableNames', {'DecileA', 'DecileB', 'Lower limit', 'A-B', 'UpperLimit', 'Pvalue'})
writetable(results, 'anovaMultipleComparison.xlsx');

%%

[p,kwtab,stats] = kruskalwallis(T.avgWithinCorrelation, T.decile)
kwResults = multcompare(stats)

results = array2table(kwResults, 'variableNames', {'DecileA', 'DecileB', 'Lower limit', 'A-B', 'UpperLimit', 'Pvalue'})
writetable(results, 'KruskalWallisMultipleComparison.xlsx');


%% try two sample t test

fs = data.avgWithinCorrelation(data.bestAlignedCorrelation <= 0.3);
fi = data.aqvgWithinCorrelation(data.bestAlignedCorrelation > 0.3);

[h,p,ci,stats] = ttest2(fi,fs, 'tail', 'right')