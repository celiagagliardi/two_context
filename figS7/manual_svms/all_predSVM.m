clear all;
data_path = '/Users/celia/Documents/two_context_data';
load(fullfile(data_path, 'tetrodes_20220316_100902_FINAL/analysis_input.mat'));
load(fullfile(data_path, 'tetrodes_20220316_100902_FINAL/analysis_results.mat'));

output_path = '/Users/celia/Documents/two_context_data/figures/SVMpred_celltype';
if ~isfolder(output_path)
    mkdir(output_path);
end


mapsData = reformatTbl(analysisInput.MapsData);
baData = reformatTbl(analysisResults.BestAligned);

fi_ba = baData(baData.isStable ==1,:);
fs_ba = baData(baData.isStable ==0,:);

fid = fopen(fullfile(output_path, 'context_pred_per_animal.txt'), 'w');
%fprintf(fid, 'Context Prediction using rates, broken down by cell type, per animal\n\n');

northDigs = {'Corr'};
southDigs = {'Geo'};
mincell = -inf;
mintrial = 4;
figure;
%% ALL cells
celltype = baData;
allCpredFi= cell(1,3); %context pred
allHpredFi = cell(1,3); % heading pred
cPredH = cell(length(7),3);
fprintf('all cells\n\n');
for d = 1:3
    dtbl = celltype(celltype.dayUsed == d,:);
    sessions = unique([dtbl(:,1), dtbl(:,end)], 'rows');
    animals = unique(sessions.animalName);
    contextPred = cell(length(animals),1);
    headingPred = cell(length(animals),1);
    for a = 1:length(animals)
        cellnames = celltype.animalSessionCellName(ismember(celltype.animalName, animals{a}) & celltype.dayUsed == d);
        stbl = mapsData(ismember(mapsData.animalName, animals{a}) & mapsData.dayUsed == d,:);
        trialTbl = unique([stbl(:,6), stbl(:,7), stbl(:,9)], 'rows'); 
        tmprv = nan(max(trialTbl.trialId), length(cellnames));
        for c = 1:length(cellnames)
            ctbl = mapsData(ismember(mapsData.animalSessionCellName, cellnames{c}),:);
            cellTrials = ctbl.trialId;
            tmprv(cellTrials, c) = ctbl.mfr;
        end
        notallnans = ~all(isnan(tmprv),1);
        tmprv = tmprv(:, notallnans);
        tmprv(isnan(tmprv)) = 0;
        validTrials = trialTbl.trialId(ismember(trialTbl.contextId, [1,2]) ...
            & ismember(trialTbl.dig, [northDigs, southDigs]));
        rv = tmprv(validTrials,:);
        if size(rv, 2) < mincell
            fprintf('animal %s day %d does not have enough cells, skipping\n', animals{a}, d);
            continue
        end

        contextClass = categorical(trialTbl.contextId(ismember(trialTbl.trialId, validTrials)));
        digClass = categorical(trialTbl.dig(ismember(trialTbl.trialId, validTrials)));
        if sum(ismember(contextClass,'1')) < 3 || sum(ismember(contextClass,'2')) < 3
            fprintf('animal %s day %d does not have enough cellular data to predict context\n', animals{a}, d);
            continue
        end

%         if sum(ismember(digClass, northDigs)) < 3 || sum(ismember(digClass, southDigs)) < 3
%             fprintf('Animal %s day %d does not have enough digs\n', animals{a}, d);
%             continue
%         end

        % encode categorical classes into numerical labels
        [~,~,Yc] = unique(contextClass);
        [~,~,Yd] = unique(digClass);

        % set up empirical prior
        classPriorC = histcounts(Yc, 1:max(Yc)+1) / numel(Yc);
        classWeightsC = 1 ./ classPriorC;
        classPriorH = histcounts(Yd, 1:max(Yd) + 1) / numel(Yd);
        classWeightsH = 1 ./ classPriorH;

       
        % initialize variables for evaluation metrics
        numSamplesC = length(Yc); %numSamplesD = length(Yd);
        confMatC = zeros(2,2);
        confMatD = zeros(2,2);

        % perform leave-one-out crossvalidation
        for iTrial = 1:numSamplesC
            %define training and testing indices
            trainIdx = [1:iTrial-1, iTrial+1:numSamplesC];
            testIdx = iTrial;

            % split data into training and test sets
            xTrain = rv(trainIdx,:);
            yTrainC = Yc(trainIdx,:);
            yTrainD = Yd(trainIdx,:);
            xTest = rv(testIdx,:);
            yTestC = Yc(testIdx,:);
            yTestD = Yd(testIdx,:);

            % Train the SVM model
            mdlC = fitcsvm(xTrain, yTrainC, 'KernelFunction', 'Linear', 'ClassNames', [1 2], 'Prior', classPriorC);
            mdlD = fitcsvm(xTrain, yTrainD, 'KernelFunction', 'Linear', 'ClassNames', [1 2], 'Prior', classPriorH);

            %predict label for test sample
            yPredC = predict(mdlC, xTest);
            yPredD = predict(mdlD, xTest);
            
            %update confuctionmatrix
            confMatC(yTestC, yPredC) = confMatC(yTestC, yPredC) + 1;
            confMatD(yTestD, yPredD) = confMatD(yTestD, yPredD) + 1;
        end

        % Evaluate the performance using confusion matrix
        accuracyC = sum(diag(confMatC)) / sum(confMatC(:));
        contextPred{a} = accuracyC;
        fprintf('Animal %s day %d context accuracy: %.2f%%\n', animals{a}, d, accuracyC * 100);
        disp('Context Confusion Matrix:');
        disp(confMatC);

        accuracyD = sum(diag(confMatD)) / sum(confMatD(:));
        cPredH{a,d} = accuracyD;
        headingPred{a} = accuracyD;
        fprintf('Animal %s day %d heading accuracy: %.2f%%\n', animals{a}, d, accuracyD * 100);
        disp('Context Confusion Matrix:');
        disp(confMatD);
    end
    allCpredFi{d} = contextPred;
    allHpredFi{d} = headingPred;
end

%%  plot the results
avg = nan(3,1);
sem1 = nan(3,1);
sem2 = nan(3,1);
fprintf('Context Prediction using ALL cell rates\n')
%fprintf(fid, 'FI cells\n');
for d = 1:3
    fprintf('Day %d\n', d)
    fprintf(fid, 'Day %d\n', d);
    cdata = cell2mat(allCpredFi{d});
    avg(d) = mean(cdata, 'omitnan');
    sem1(d) = std(cdata, 'omitnan') ./ sqrt(sum(~isnan(cdata)));
    sem2(d) = sqrt((mean(cdata, 'omitnan')*(1-mean(cdata, 'omitnan'))) ./ sum(~isnan(cdata)));

    [h,p, ci, stats] = ttest(cdata, .5, 'tail', 'right');
    %fprintf('\tContext Prediction on day %d p val = %.3f\n', d, p);

    fprintf('\t avg = %.4f, sem = %.4f, t(%d) = %.3f, p = %.5f, using %d animals\n',  avg(d), sem1(d), stats.df, stats.tstat, p, sum(~isnan(cdata)));
    %fprintf(fid, '\t avg = %.4f, sem = %.4f, t(%d) = %.3f, p = %.5f, using %d animals\n',  avg(d), sem1(d), stats.df, stats.tstat, p, sum(~isnan(cdata)));
end

sem = sem1;

% Example data 
subplot(1,2,1);
b = bar(avg, 'grouped');
hold on
% Find the number of groups and the number of bars in each group
[ngroups, nbars] = size(avg);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, avg(:,i), sem(:,i), 'k', 'linestyle', 'none');
end


data1 = cell2mat(allCpredFi{1});
data2 = cell2mat(allCpredFi{2});
data3 = cell2mat(allCpredFi{3});

scatter(repmat(1, length(data1),1), data1, 'ro', 'xjitter', 'rand')
scatter(repmat(2, length(data2),1), data2, 'bo', 'xjitter', 'rand')
scatter(repmat(3, length(data3),1), data3, 'co', 'xjitter', 'rand')

hold off
ylim([0 1])
ylabel('Prediction Accuracy');
xlabel('Day');
title('ALL cell context prediction');

% heading prediction      

avg = nan(3,1);
sem1 = nan(3,1);
sem2 = nan(3,1);
fprintf('Heading Prediction using ALL cell rates\n')
fprintf(fid, 'FI cells\n');
for d = 1:3
    fprintf('Day %d\n', d)
    fprintf(fid, 'Day %d\n', d);
    cdata = cell2mat(allHpredFi{d});
    avg(d) = mean(cdata, 'omitnan');
    sem1(d) = std(cdata, 'omitnan') ./ sqrt(sum(~isnan(cdata)));
    sem2(d) = sqrt((mean(cdata, 'omitnan')*(1-mean(cdata, 'omitnan'))) ./ sum(~isnan(cdata)));

    %[h,p, ci, stats] = ttest(cdata, .5, 'tail', 'right');
    %fprintf('\tContext Prediction on day %d p val = %.3f\n', d, p);
    [p,h] = signrank(cdata, .5, 'tail', 'right');

    fprintf('\t avg = %.4f, sem = %.4f, h = %d, p = %.5f, using %d animals\n',  avg(d), sem1(d), h, p, sum(~isnan(cdata)));
    %fprintf(fid, '\t avg = %.4f, sem = %.4f, t(%d) = %.3f, p = %.5f, using %d animals\n',  avg(d), sem1(d), stats.df, stats.tstat, p, sum(~isnan(cdata)));
end

sem = sem1;

% Example data 
subplot(1,2,2);
b = bar(avg, 'grouped');
hold on
% Find the number of groups and the number of bars in each group
[ngroups, nbars] = size(avg);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, avg(:,i), sem(:,i), 'k', 'linestyle', 'none');
end


data1 = cell2mat(allHpredFi{1});
data2 = cell2mat(allHpredFi{2});
data3 = cell2mat(allHpredFi{3});

scatter(repmat(1, length(data1),1), data1, 'ro', 'xjitter', 'rand')
scatter(repmat(2, length(data2),1), data2, 'bo', 'xjitter', 'rand')
scatter(repmat(3, length(data3),1), data3, 'co', 'xjitter', 'rand')

hold off
ylim([0 1])
ylabel('Prediction Accuracy');
xlabel('Day');
title('ALL cell heading prediction');

sgtitle('with empirical prior');
       