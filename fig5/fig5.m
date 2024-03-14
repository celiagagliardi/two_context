
%% load data
% run compute_coherency_table and then summarize_coherency_table before you
% run the rest of this script. 
% 
% filepath = fileparts(mfilename('fullpath'));
% cd(filepath);
% compute_coherency_table;
% summarize_coherency_table;
% 
%% Figure 5C
data = readtable('output/coherency_table.csv');
summary = readtable('output/coherency_summary.csv');
coh = summary(summary.FIisCoherent == 1,:);

trialPair = coh(randi(size(coh,1),1,1),:);

t = data(ismember(data.animal, trialPair.animal) & data.day == trialPair.day & data.trialA == trialPair.trialA & data.trialB == trialPair.trialB,:);

bins = 0:.01:1;
histogram(summary.FI_rot1, bins); hold on;
histogram(summary.FI_rot2,bins);
histogram(summary.FI_rot3,bins);
histogram(summary.FI_rot4, bins);
xline(.25)
xlabel('Percent Cells'); % multiply by 100 in illustrator
ylabel('Trial Pairs');

%% Figure 5D

animals = unique(summary.animal);

fi_within_trialPairs = cell(length(animals),3);
fi_across_trialPairs = cell(length(animals),3);
fs_within_trialPairs = cell(length(animals),3);
fs_across_trialPairs = cell(length(animals),3);

fi_within_rotdist1 = cell(length(animals),3);
fi_across_rotdist1 = cell(length(animals),3);
fs_within_rotdist1 = cell(length(animals),3);
fs_across_rotdist1 = cell(length(animals),3);

fi_within_rotdist0 = cell(length(animals),3);
fi_across_rotdist0 = cell(length(animals),3);
fs_within_rotdist0 = cell(length(animals),3);
fs_across_rotdist0 = cell(length(animals),3);

fi_within_rotdistALL = cell(length(animals),3);
fi_across_rotdistALL = cell(length(animals),3);
fs_within_rotdistALL = cell(length(animals),3);
fs_across_rotdistALL = cell(length(animals),3);

fi_within_rotdist = cell(length(animals),3); % per trial pair
fi_across_rotdist = cell(length(animals),3);
fs_within_rotdist = cell(length(animals),3);
fs_across_rotdist = cell(length(animals),3);

for a = 1:length(animals)
    for d = 1:3
        dtbl = summary(ismember(summary.animal, animals{a}) & summary.day == d,:);
        within = dtbl(ismember(dtbl.comparison, 'within'),:);
        fi_within_trialPairs{a,d} = mean(within.FIisCoherent, 'omitnan');
        fs_within_trialPairs{a,d} = mean(within.FSisCoherent, 'omitnan');
        
        fi_within_rotdist1{a,d} = mean(table2array(within(within.FIisCoherent == 1, 9:12)),1, 'omitnan');
        fs_within_rotdist1{a,d} = mean(table2array(within(within.FSisCoherent == 1, 14:17)),1, 'omitnan');
        fi_within_rotdist0{a,d} = mean(table2array(within(within.FIisCoherent == 0, 9:12)),1, 'omitnan');
        fs_within_rotdist0{a,d} = mean(table2array(within(within.FSisCoherent == 0, 14:17)),1, 'omitnan');

        fi_within_rotdistALL{a,d} = mean(table2array(within(:, 9:12)), 1, 'omitnan');
        fs_within_rotdistALL{a,d} = mean(table2array(within(:, 14:17)),1, 'omitnan');

        fi_within_rotdist{a,d} = table2array(within(:, 9:12));
        fs_within_rotdist{a,d} = table2array(within(:, 14:17));


        across = dtbl(ismember(dtbl.comparison, 'across'),:);
        fi_across_trialPairs{a,d} = mean(across.FIisCoherent, 'omitnan');
        fs_across_trialPairs{a,d} = mean(across.FSisCoherent, 'omitnan');
        
        fi_across_rotdist1{a,d} = mean(table2array(across(across.FIisCoherent == 1, 9:12)),1, 'omitnan');
        fs_across_rotdist1{a,d} = mean(table2array(across(across.FSisCoherent == 1, 14:17)),1, 'omitnan');
        
        fi_across_rotdist0{a,d} = mean(table2array(across(across.FIisCoherent == 0, 9:12)),1, 'omitnan');
        fs_across_rotdist0{a,d} = mean(table2array(across(across.FSisCoherent == 0, 14:17)),1, 'omitnan');
        
        fi_across_rotdistALL{a,d} = mean(table2array(across(:, 9:12)),1, 'omitnan');
        fs_across_rotdistALL{a,d} = mean(table2array(across(:, 14:17)),1, 'omitnan');

        fi_across_rotdist{a,d} = table2array(across(:, 9:12));
        fs_across_rotdist{a,d} = table2array(across(:, 14:17));

    end
end



fprintf('Rotation distribution of all coherent trial pairs per trial pair')
figure
set(gcf, 'position', [1, 467, 1728, 510]);
for d = 1:3

    avg = nan(4,4);
    sem = nan(4,4);

    fprintf('Day %d\n', d)
    %fprintf(fid, '\tDay %d\n', d);

    cdata = cell2mat(fi_within_rotdist(:,d));
    avg(1,:) = mean(cdata, 'omitnan');
    sem(1,:) = std(cdata, 'omitnan') ./ sqrt(sum(~isnan(cdata)));

    cdata = cell2mat(fi_across_rotdist(:,d));
    avg(2,:) = mean(cdata, 'omitnan');
    sem(2,:) = std(cdata, 'omitnan') ./ sqrt(sum(~isnan(cdata)));
     
    cdata = cell2mat(fs_within_rotdist(:,d));
    avg(3,:) = mean(cdata, 'omitnan');
    sem(3,:) = std(cdata, 'omitnan') ./ sqrt(sum(~isnan(cdata)));

    cdata = cell2mat(fs_across_rotdist(:,d));
    avg(4,:) = mean(cdata, 'omitnan');
    sem(4,:) = std(cdata, 'omitnan') ./ sqrt(sum(~isnan(cdata)));

    ax(d) = subplot(1,3,d);
    b = bar(avg, 'grouped');
    hold on
    [ngroups, nbars] = size(avg);
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    for i = 1:nbars
        % Calculate center of each bar
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x, avg(:,i), sem(:,i), 'k', 'linestyle', 'none');
    end
    yline(.25, 'r-')
    xticklabels({'FI Within', 'FI Across', 'FS Within', 'FS Across'});
end
sgtitle('Rotation distribution in all trial pairs per trial pair')
linkaxes(ax, 'xy');

%% Figure 5E
 
fiwithin = cell(length(animals),3);
fiacross = cell(length(animals),3);
fswithin = cell(length(animals),3);
fsacross = cell(length(animals),3);
for a = 1:length(animals)
    for d = 1:3
        dtbl = summary(ismember(summary.animal, animals{a}) & summary.day == d,:);
        fiwithin{a,d} = [dtbl.FI_rot1(ismember(dtbl.comparison, 'within')), ...
            dtbl.FI_rot2(ismember(dtbl.comparison, 'within'))];
        fswithin{a,d} = [dtbl.FS_rot1(ismember(dtbl.comparison, 'within')), ...
            dtbl.FS_rot2(ismember(dtbl.comparison, 'within'))];

        fiacross{a,d} = [dtbl.FI_rot1(ismember(dtbl.comparison, 'across')), ...
            dtbl.FI_rot2(ismember(dtbl.comparison, 'across'))];
        fsacross{a,d} = [dtbl.FS_rot1(ismember(dtbl.comparison, 'across')), ...
            dtbl.FS_rot2(ismember(dtbl.comparison, 'across'))];
    end
end

figure;
ax(1) = subplot(2,2,1);

data1 = cell2mat(fiwithin(:,1));
data2 = cell2mat(fiwithin(:,2));
data3 = cell2mat(fiwithin(:,3));
avg1= mean(data1);
avg2 = mean(data2);
avg3 = mean(data3);

sem1 = std(data1) ./ sqrt(size(data1,1));
sem2 = std(data2) ./ sqrt(size(data2,1));
sem3 = std(data3) ./ sqrt(size(data3,1));

avg = [avg1; avg2; avg3];
sem = [sem1; sem2; sem3];

b = bar(avg, 'grouped', 'basevalue', .25);
hold on
[ngroups, nbars] = size(avg);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, avg(:,i), sem(:,i), 'k', 'linestyle', 'none');
end
title('FI Within');
xticklabels({'Day 1', 'Day 2', 'Day 3'})
ylabel('Cell Proportions');

ax(2) = subplot(2,2,2);
data1 = cell2mat(fiacross(:,1));
data2 = cell2mat(fiacross(:,2));
data3 = cell2mat(fiacross(:,3));
avg1= mean(data1);
avg2 = mean(data2);
avg3 = mean(data3);

sem1 = std(data1) ./ sqrt(size(data1,1));
sem2 = std(data2) ./ sqrt(size(data2,1));
sem3 = std(data3) ./ sqrt(size(data3,1));

avg = [avg1; avg2; avg3];
sem = [sem1; sem2; sem3];

b = bar(avg, 'grouped', 'basevalue', .25);
hold on
[ngroups, nbars] = size(avg);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, avg(:,i), sem(:,i), 'k', 'linestyle', 'none');
end
title('FI Across');
xticklabels({'Day 1', 'Day 2', 'Day 3'})
ylabel('Cell Proportions');

ax(3) = subplot(2,2,3);

data1 = cell2mat(fswithin(:,1));
data2 = cell2mat(fswithin(:,2));
data3 = cell2mat(fswithin(:,3));
avg1= mean(data1, 'omitnan');
avg2 = mean(data2, 'omitnan');
avg3 = mean(data3, 'omitnan');

sem1 = std(data1, 'omitnan') ./ sqrt(sum(~isnan(data1)));
sem2 = std(data2, 'omitnan') ./ sqrt(sum(~isnan(data2)));
sem3 = std(data3, 'omitnan') ./ sqrt(sum(~isnan(data3)));

avg = [avg1; avg2; avg3];
sem = [sem1; sem2; sem3];

b = bar(avg, 'grouped', 'basevalue', .25);
hold on
[ngroups, nbars] = size(avg);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, avg(:,i), sem(:,i), 'k', 'linestyle', 'none');
end
title('FS Within');
xticklabels({'Day 1', 'Day 2', 'Day 3'})
ylabel('Cell Proportions');

ax(4) = subplot(2,2,4);
data1 = cell2mat(fsacross(:,1));
data2 = cell2mat(fsacross(:,2));
data3 = cell2mat(fsacross(:,3));
avg1= mean(data1, 'omitnan');
avg2 = mean(data2, 'omitnan');
avg3 = mean(data3, 'omitnan');

sem1 = std(data1, 'omitnan') ./ sqrt(sum(~isnan(data1)));
sem2 = std(data2, 'omitnan') ./ sqrt(sum(~isnan(data2)));
sem3 = std(data3, 'omitnan') ./ sqrt(sum(~isnan(data3)));

avg = [avg1; avg2; avg3];
sem = [sem1; sem2; sem3];

b = bar(avg, 'grouped', 'basevalue', .25);
hold on
[ngroups, nbars] = size(avg);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, avg(:,i), sem(:,i), 'k', 'linestyle', 'none');
end
title('FS Across');
xticklabels({'Day 1', 'Day 2', 'Day 3'})
ylabel('Cell Proportions');
legend('1st BMR', '2nd BMR')

linkaxes(ax, 'xy')

%% Figure 5F

data = readtable('output/coherency_table.csv');
animals = unique(data.animal);

fi_within = cell(length(animals),3);
fi_across = cell(length(animals),3);
fs_within = cell(length(animals),3);
fs_across = cell(length(animals),3);


for a = 1:length(animals)
    for d = 1:3
        dtbl = data(ismember(data.animal, animals{a}) & data.day == d,:);
        pwTrials = unique([dtbl.trialA, dtbl.trialB], 'rows');
        contextChanges = [mod(pwTrials(:,1),2), mod(pwTrials(:,2),2)];
        withinContext = contextChanges(:,1)== contextChanges(:,2);

        fi_rotidx = nan(size(pwTrials,1),4);
        fs_rotidx = nan(size(pwTrials,1),4);

        for pw = 1:size(pwTrials,1)
            pwTbl = dtbl(dtbl.trialA == pwTrials(pw,1) & dtbl.trialB == pwTrials(pw,2),:);

            % FI CELLS

            bmr = pwTbl.bmr(~isnan(pwTbl.bmr) & pwTbl.isStable == 1);

            if length(bmr) > 4

                n = histcounts(bmr, 1:5, 'normalization', 'probability');
                [~, idx_n] = sort(n, 'descend') ;
                fi_rotidx(pw,:) = idx_n;
            else
                fprintf('%s day %d trial pair %d and %d do not have enough FI cells, skipping\n', animals{a}, d, pwTrials(pw,1), pwTrials(pw,2));
            end


            % FS CELLS
            bmr = pwTbl.bmr(~isnan(pwTbl.bmr) & pwTbl.isStable == 0);
            if length(bmr) > 4
                n = histcounts(bmr, 1:5, 'normalization', 'probability');
                [~, idx_n] = sort(n, 'descend');
                fs_rotidx(pw,:) = idx_n;
            else
                fprintf('%s day %d trial pair %d and %d do not have enough FS cells, skipping\n', animals{a}, d, pwTrials(pw,1), pwTrials(pw,2));
            end


        end
        fi_within{a,d} = fi_rotidx(withinContext,:);
        fi_across{a,d} = fi_rotidx(~withinContext,:);

        fs_within{a,d} = fs_rotidx(withinContext,:);
        fs_across{a,d} = fs_rotidx(~withinContext,:);

    end
end
% plot

fiw_all = cell(length(animals), 3);
fia_all = fiw_all;
fsw_all = fiw_all;
fsa_all = fiw_all;


withinContext = data.contextA == data.contextB;


for d = 1:3
    
    fiwcells = sum(~isnan(data.bmr(data.day == d & data.isStable == 1 & data.contextA == data.contextB)));
    fiacells = sum(~isnan(data.bmr(data.day == d & data.isStable == 1 & data.contextA ~= data.contextB)));

    fswcells = sum(~isnan(data.bmr(data.day == d & data.isStable == 0 & data.contextA == data.contextB)));
    fsacells = sum(~isnan(data.bmr(data.day == d & data.isStable == 0 & data.contextA ~= data.contextB)));


    fiw = nan(4,4, length(animals));
    fia = fiw;
    fsw = fiw;
    fsa = fiw;   
       
    for a = 1:length(animals)
        cdata = cell2mat(fi_within(a,d));
        fiwn = sum(~isnan(data.bmr(data.day == d & data.isStable == 1 & data.contextA == data.contextB & ismember(data.animal, animals{a}))));

        if ~all(isnan(cdata), 'all')
            tmp = zeros(4,4);
            for c = 1:size(cdata,1)
                tmp(cdata(c,1), cdata(c,2)) = tmp(cdata(c,1), cdata(c,2)) + 1;
            end
            fiw(:,:,a) = tmp ./ sum(tmp, 'all');
            fiw_all{a,d} = fiw(:,:,a);
        end

        cdata = cell2mat(fi_across(a,d));
        fian = sum(~isnan(data.bmr(data.day == d & data.isStable == 1 & data.contextA ~= data.contextB & ismember(data.animal, animals{a}))));


        if ~all(isnan(cdata), 'all')
            tmp = zeros(4,4);
            for c = 1:size(cdata,1)
                tmp(cdata(c,1), cdata(c,2)) = tmp(cdata(c,1), cdata(c,2)) + 1;
            end
            fia(:,:,a) = tmp ./ sum(tmp, 'all');
            fia_all{a,d} = fia(:,:,a);
        end

        cdata = cell2mat(fs_within(a,d));
        fswn = sum(~isnan(data.bmr(data.day == d & data.isStable == 0 & data.contextA == data.contextB & ismember(data.animal, animals{a}))));

        if ~all(isnan(cdata), 'all')

            tmp = zeros(4,4);
            for c = 1:size(cdata,1)
                tmp(cdata(c,1), cdata(c,2)) = tmp(cdata(c,1), cdata(c,2)) + 1;
            end
            fsw(:,:,a) = tmp ./ sum(tmp, 'all');
            fsw_all{a,d} = fsw(:,:,a);
        end

        cdata = cell2mat(fs_across(a,d));
        fsan = sum(~isnan(data.bmr(data.day == d & data.isStable == 0 & data.contextA ~= data.contextB & ismember(data.animal, animals{a}))));

        if ~all(isnan(cdata), 'all')

            tmp = zeros(4,4);
            for c = 1:size(cdata,1)
                tmp(cdata(c,1), cdata(c,2)) = tmp(cdata(c,1), cdata(c,2)) + 1;
            end
            fsa(:,:,a) = tmp ./ sum(tmp, 'all');
            fsa_all{a,d} = fsa(:,:,a);
        end

    end

    fiwavg = mean(fiw, 3, 'omitnan');
    fiaavg = mean(fia, 3, 'omitnan');
    fswavg = mean(fsw, 3, 'omitnan');
    fsaavg = mean(fsa, 3, 'omitnan');

    % omit diagonal as per reviewer 3 request
    fiwavg(logical(eye(4))) = nan;
    fiaavg(logical(eye(4))) = nan;
    fswavg(logical(eye(4))) = nan;
    fsaavg(logical(eye(4))) = nan;


    figure
    subplot(2,2,1)
    hfig = imagesc(fiwavg);
    set(hfig, 'AlphaData', ~isnan(fiwavg))
    colorbar;
    ylabel('First BMR');
    xlabel('Second BMR');
    xticks([1 2 3 4]);
    yticks([1 2 3 4]);
    yticklabels({'0', '90', '180', '270'});
    xticklabels({'0', '90', '180', '270'});
    title('fi within');

    subplot(2,2,2)
    hfig = imagesc(fiaavg);
    set(hfig, 'AlphaData', ~isnan(fiaavg))
    colorbar;
    ylabel('First BMR');
    xlabel('Second BMR');
    xticks([1 2 3 4]);
    yticks([1 2 3 4]);
    yticklabels({'0', '90', '180', '270'});
    xticklabels({'0', '90', '180', '270'});
    title('fi across');

    subplot(2,2,3)
    hfig = imagesc(fswavg);
    set(hfig, 'AlphaData', ~isnan(fswavg))
    colorbar;
    ylabel('First BMR');
    xlabel('Second BMR');
    xticks([1 2 3 4]);
    yticks([1 2 3 4]);
    yticklabels({'0', '90', '180', '270'});
    xticklabels({'0', '90', '180', '270'});
    title('fs within');

    subplot(2,2,4)
    hfig = imagesc(fsaavg);
    set(hfig, 'AlphaData', ~isnan(fsaavg))
    colorbar;
    ylabel('First BMR');
    xlabel('Second BMR');
    xticks([1 2 3 4]);
    yticks([1 2 3 4]);
    yticklabels({'0', '90', '180', '270'});
    xticklabels({'0', '90', '180', '270'});
    title('fs across');

    sgtitle(sprintf('day %d', d));

end

fiw_all(cellfun(@isempty, fiw_all)) = {nan(4,4)};
fia_all(cellfun(@isempty, fia_all)) = {nan(4,4)};
fsw_all(cellfun(@isempty, fsw_all)) = {nan(4,4)};
fsa_all(cellfun(@isempty, fsa_all)) = {nan(4,4)};

%tabulate

for d = 1:3
    figure
    cdata = cell2mat(fi_within(:,d));
    fprintf('fi within\n')
    t = array2table(tabulate(cdata(:,1)), 'VariableNames', {'Value', 'Count', 'Frequency'});
    ax(1) = subplot(2,2,1);
    bar(t.Value, t.Frequency)
    title('fi within');


    cdata = cell2mat(fi_across(:,d));
    fprintf('fi across\n')
    t = array2table(tabulate(cdata(:,1)), 'VariableNames', {'Value', 'Count', 'Frequency'});
    ax(2) = subplot(2,2,2);
    bar(t.Value, t.Frequency)
    title('fi across');

    cdata = cell2mat(fs_within(:,d));
    fprintf('fs within\n')
    t = array2table(tabulate(cdata(:,1)), 'VariableNames', {'Value', 'Count', 'Frequency'});
    ax(3) = subplot(2,2,3);
    bar(t.Value, t.Frequency)
    title('fs within');

    cdata = cell2mat(fs_across(:,d));
    fprintf('fs across\n')
    t = array2table(tabulate(cdata(:,1)), 'VariableNames', {'Value', 'Count', 'Frequency'});
    ax(4) = subplot(2,2,4);
    bar(t.Value, t.Frequency)
    title('fs across');

    sgtitle(sprintf('Day %d 1st BMR Marginal Probabilities', d));
    linkaxes(ax, 'xy');
end

