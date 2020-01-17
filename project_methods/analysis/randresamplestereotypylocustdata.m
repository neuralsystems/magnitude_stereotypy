function randresamplestereotypylocustdata()
if ~exist('data\locust_bln1\bln1_stereotypy_resample.mat', 'file')
    rng('default');
    load('data\locust_bln1\bln1_stereotypy.mat')
    nTrial = 1e5;
    stereotypy = cell(1, nTrial);
    for iTrial = 1:nTrial
        data = reshape(meanfiringNorm(randperm(numel(meanfiringNorm))), 6, []);
        stereotypy{iTrial} = computeindividualpairstereotypy(data);
    end
    meanStereo = cellfun(@computeultimatemean, stereotypy);
    % random resample < mean actual
    % REPORTED VALUE %
    p = 1 - sum(meanStereo < mean(all_values)) / nTrial;
    % left-tailed ttest, random resample vs mean actual
    [~, p1] = ttest(meanStereo, mean(all_values), 'tail', 'left');
    % two-sided ttest, random resample vs 0
    [~, p2] = ttest(meanStereo, 0);
    % two-sided ttest, random resample vs mean actual
    [~, p3] = ttest(meanStereo, mean(all_values));
    save('data\locust_bln1\bln1_stereotypy_resample.mat', 'p*', 'meanStereo', 'nTrial');
else
    load('data\locust_bln1\bln1_stereotypy_resample.mat');
end
fprintf('Reported p-value (random resample < mean actual) = %2.4g\n', p);
fprintf('p-value (left-tailed ttest, random resample vs mean actual) = %2.4g\n', p1);
fprintf('p-value (two-sided ttest, random resample vs 0) = %2.4g\n', p2);
fprintf('p-value (two-sided ttest, random resample vs mean actual) = %2.4g\n', p3);
end