clear all;

data_path = '/Users/celia/Documents/two_context_data';
load(fullfile(data_path, "ratemaps-2021_12_1_12_34_49.mat"));

output_path = '/Users/celia/Documents/two_context_data/figures/trajectories';
if ~isfolder(output_path)
    mkdir(output_path);
end

% plot the trajectory of the animal and pick one

for a = 1:length(ratemaps.data)
    for d = 1:length(ratemaps.data(a).session)
        figure;
        for t = 1:length(ratemaps.data(a).session(d).trial)
            x = ratemaps.data(a).session(d).trial(t).x;
            y = ratemaps.data(a).session(d).trial(t).y;
            subplot(2,6,t);
            plot(x,y); axis ij; daspect([1 1 1]);
            title(sprintf('Trial %d', ratemaps.data(a).session(d).trial(t).trialNum));
        end
        sgtitle({sprintf('%s', ratemaps.data(a).animal), sprintf('Day %d', d)}, ...
            'interpreter', 'none');
        saveas(gcf, fullfile(output_path, sprintf('%s_%s_trajectories.pdf', ratemaps.data(a).animal, ratemaps.data(a).session(d).name)))
    end
end

close all;


    