%% Figure S13B
data = readtable("figure_S13B_alternate_cellreg_count_results.xlsx");
panels = [{'Day 1'}, {'Day 2'}; {'Day 1'} {'Day 3'}; {'Day 2'}, {'Day 3'}];

for p = 1:size(panels,1);
    fprintf('Comparing %s with %s\n', panels{p,1}, panels{p,2});
    dtbl = data(ismember(data.groupLabelA, panels(p,1)) & ismember(data.groupLabelB, panels(p,2)),:);
    X = abs(dtbl.bestCorrelationChange);
    G = dtbl.binIdA;

    figure
    [pval, anovatab, stats] = kruskalwallis(X,G);
    ylabel('Absolute correlation difference');
    xlabel('Decile');
    figure
    T = multcompare(stats);
    title('Multiple comparison of mean ranks');
    xlabel('Mean Ranks');
    ylabel('Decile');
    tbl = array2table(T, 'VariableNames', {'DecileA', 'DecileB', 'LowerLimit', 'A-B', 'Upper Limit', 'P-value'});
    %writetable(tbl, sprintf('%s_to_%s_multiplecomparison.xlsx', panels{p,1}, panels{p,2}));

    title(sprintf('%s to %s', panels{p,1}, panels{p,2}))
end

%% Try two way ranked anova.
corrdata = abs(data.bestCorrelationChange);
[~, corrrank] = sort(corrdata, 'ascend');
rankData = [data, table(corrrank, 'VariableName', {'RankedData'})];
fprintf('Ranked Anova\n');

X = rankData.RankedData;
g1 = rankData.binIdA;
% g2 = nan(size(rankData,1),1);
% g2(ismember(rankData.groupLabelA))si
g2 = rankData.groupLabelA;
g3 = rankData.groupLabelB;

[p,tbl,stats,terms] = anovan(X, {g1, g2, g3}, 'model', 'interaction', 'varnames', {'g1', 'g2', 'g3'});

g2 =  strcat(rankData.groupLabelA, rankData.groupLabelB);

[p,tbl,stats,terms] = anovan(X, {g1, g2}, 'model', 'interaction', 'varnames', {'g1', 'g2'});
[results,~,~,gnames] = multcompare(stats,"Dimension",[1 2]);
tbl = array2table(results,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tbl.("Group A")=gnames(tbl.("Group A"));
tbl.("Group B")=gnames(tbl.("Group B"))
