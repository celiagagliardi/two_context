data_path = '/Users/celia/Documents/two_context_data';
calin = load(fullfile(data_path, 'calcium_20220516_170618/analysis_input.mat'));
maps = reformatTbl(calin.analysisInput.MapsData);
calout = load(fullfile(data_path, 'calcium_20220516_170618/analysis_results.mat'));
ba = reformatTbl(calout.analysisResults.BestAligned);
ficells = ba.animalSessionCellName(ba.isStable == 1);
tbl = maps(ismember(maps.animalSessionCellName, ficells),:);

panels = [1,2;1,3; 2, 3]; % col 1 is train day and col2 is predict day
northDigs = {'Corr', 'Feat'}; southDigs = {'Geo', 'Wrong'};

%% 
X = cell(1,3);
for t = 1:size(panels,1)
    fprintf('Train Day %d, predict day %d\n', panels(t,1), panels(t,2));
    cellnames = intersect(tbl.registeredCellName(tbl.dayUsed == panels(t,1)), tbl.registeredCellName(tbl.dayUsed == panels(t,2)));
    
    cPred = nan(length(cellnames),1);
    for c = 1:length(cellnames)
        % train
        trainData = tbl(ismember(tbl.registeredCellName, cellnames{c}) & tbl.dayUsed == panels(t,1),:);
        validData = trainData(ismember(trainData.dig, [northDigs, southDigs]), :);
        if isempty(validData)
            fprintf('not enough training data, skipping\n');
            continue
        end
        digClass = categorical(validData.dig);
        digClass(ismember(digClass, northDigs)) = 'North';
        digClass(ismember(digClass, southDigs)) = 'South';
        rv = cellfun(@centerOutRadians, validData.map);
        trainVector = [sin(rv), cos(rv)];
        mdl = fitcsvm(trainVector, digClass, 'standardize', true, 'kernelfunction', 'polynomial'); %% NO CROSS VALIDATION
       
        % predict
        testData = tbl(ismember(tbl.registeredCellName, cellnames{c}) & tbl.dayUsed == panels(t,2),:);
        validData = testData(ismember(testData.dig, [northDigs, southDigs]), :);
        testDig = categorical(validData.dig);
        testDig(ismember(testDig, northDigs)) = 'North';
        testDig(ismember(testDig, southDigs)) = 'South';
        rv = cellfun(@centerOutRadians, validData.map);
        testVector = [sin(rv), cos(rv)];
      
        label = predict(mdl, testVector);
        pred_acc = mean(testDig == label);
        cPred(c) = pred_acc;
    end

    
X{t} = cPred;
    avg = mean(cPred);
    sem = std(cPred) ./ sqrt(length(cPred));
    fprintf('\t avg = %.3f, sem = %.3f\n', avg, sem);
    p = signrank(cPred, .5, 'tail', 'right');
    fprintf('\t\t sign rank: p =  %.3f\n', p);
     [~,p,ci, stats] = ttest(cPred, .5, 'tail', 'right');
        fprintf('One sample t test: t(%d) = %.3f, p = %.3f\n', stats.df, stats.tstat, p);
     h = kstest(cPred);
    if h
        fprintf('\t\tdata not normal, recommend sign rank test\n');
    end

end

X = X';
l = cellfun(@length, X);
g = [repmat(1, l(1),1); repmat(2,l(2),1); repmat(3,l(3),1)];
figure;
boxplot(cell2mat(X), g);
yline(.5);
ylim([0 1.2]);
title('Per cell')


