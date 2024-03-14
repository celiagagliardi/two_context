data_path = '/Users/celia/Documents/two_context_data'; % change this

load(fullfile(data_path, 'ratemaps-2021_12_1_12_34_49_K1_d4.mat'))
minTrial = 4;
animals = length(extractfield(ratemaps.data, 'animal'));
within = cell(length(animals),3);
across = cell(length(animals),3);
rmDay = cell(length(animals), 3);

for a = 1:length(ratemaps.data)
    switch ratemaps.data(a).animal
        case 'AK42_CA1'
            sessionName = {'d7', 'd8', 'd9'};
        case 'AK74_CA1'
            sessionName = {'d1', 'd2', 'd3'};
        case 'CMG159_recut'
            sessionName = {'s6', 's7', 's8'};
        case 'HGY1_recut'
            sessionName = {'d1'};
        case 'JJ9_CA1'
            sessionName = {'d1', 'd2', 'd3'};
        case 'K1_CA1'
            sessionName = {'d1', 'd2', 'd4'};
        case 'MG1'
            sessionName = {'d6', 'd7', 'd11'};
    end
    
    sess = find(ismember(extractfield(ratemaps.data(a).session, 'name'), sessionName));
    
    for s = 1:length(sess)
        trials = extractfield(ratemaps.data(a).session(sess(s)).trial, 'trialNum');
        contexts = extractfield(ratemaps.data(a).session(sess(s)).trial, 'context');
        trList = [find(contexts == 1), find(contexts == 2)];
        digs = extractfield(ratemaps.data(a).session(sess(s)).trial, 'dig');
        allcells = findUniqueCells(ratemaps.data(a).session(sess(s)));
        withinCell = nan(length(allcells),1);
        acrossCell = nan(length(allcells),1);
        rmCell = nan(12,12,length(allcells));
        for c = 1:length(allcells)
            mfr = calculateMfrVector(ratemaps.data(a).session(sess(s)), allcells{c});
            trialsAvail = trials(ismember(digs, {'Corr', 'Geo', 'Feat', 'Wrong'}) & ~isnan(mfr));
            if length(trialsAvail) < minTrial
                continue;
            end
            
            com = nchoosek(trialsAvail,2);
            totalWithin = nan(size(com,1),1);
            totalAcross = nan(size(com,1),1);
            for p = 1:size(com,1)
               
                mfrdiff = abs(mfr(com(p,1)) - mfr(com(p,2)));

                indA = find(trList == com(p,1));
                indB = find(trList == com(p,2));
                rmCell(indA, indB, c) = mfrdiff;
                rmCell(indB, indA, c) = mfrdiff;
                    
                if contexts(com(p,1)) == contexts(com(p,2))
                    totalWithin(p) = mfrdiff;
                else
                    totalAcross(p) = mfrdiff;
                end
            end
            
            if all(isnan(totalWithin)) || all(isnan(totalAcross))
                continue;
            end
            
            withinCell(c) = nanmean(totalWithin);
            acrossCell(c) = nanmean(totalAcross);
            
        end
        
        within{a,s} = withinCell;
        across{a,s} = acrossCell;
        rmDay{a,s} = rmCell;
    end
    
end
clim = [.096, .337];
for d = 1:3
    dayRm = rmDay(:,d);
    allRm = cat(3,dayRm{:});
    avgRm = nanmean(allRm,3);
    figure;
    ax = plotRRM(avgRm);
    %colorbar;
    fprintf('day %d clim is %.3f, %.3f\n', d, get(gca, 'clim'))
    title(sprintf('day %d', d));
    caxis(clim)
end

colorbar

