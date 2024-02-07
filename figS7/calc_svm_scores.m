function T = calc_svm_scores(hMdl, rv, truelabels)

testInds = nan(size(rv,1),1);
for t = 1:size(testInds)
    testInds(t) = find(test(hMdl.Partition,t));
end

xTest = rv(testInds,:);
yTest = truelabels(testInds,:);
compactMdl = hMdl.Trained{1};
[label,score] = predict(compactMdl,xTest);
T = table(yTest,label,score(:,2),'VariableNames',...
    {'TrueLabel','PredictedLabel','Score'});