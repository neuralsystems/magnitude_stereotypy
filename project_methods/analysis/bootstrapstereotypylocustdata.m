function bootstrapstereotypylocustdata()
rng('default');
load('data\bln1_stereotypy.mat')
nTrial = 1e5;
stereotypy = cell(1, nTrial);
for iTrial = 1:nTrial
    data = reshape(meanfiringNorm(randperm(numel(meanfiringNorm))), 6, []);
    stereotypy{iTrial} = computeindividualpairstereotypy(data);
end
meanStereo = cellfun(@computeultimatemean, stereotypy);
p = 1 - sum(meanStereo < mean(all_values)) / 1e5;
disp(p)
[~, p1] = ttest(meanStereo, mean(all_values), 'tail', 'left');
disp(p1)
[~, p2] = ttest(meanStereo, 0);
disp(p2)
[~, p3] = ttest(meanStereo, mean(all_values));
disp(p3)
save('data\bln1_stereotypy_boot.mat');
end