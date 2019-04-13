function run = savepnvectorcorrelationvsstereotypyplot(pathResult, nSeed)

%**********************************************************************%
% Author: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

try
    % load existing data file
    if exist([pathResult, 'raw_data.mat'], 'file')
        load([pathResult, 'raw_data.mat']);
    else
        error('raw_data.mat file doesn''t exist');
    end
    if exist([pathResult, 'spike_data.mat'], 'file')
        load([pathResult, 'spike_data.mat']);
        if ~exist('spikeData', 'var')
            error('spikeData variable not found in spike_data.mat');
        end
    else
        error('spike_data.mat file doesn''t exist');
    end
    % define layer types and titles
    pathPlot = createresultfolder([pathResult, 'correlation\'], true);
    idLayerCorr = {'vopn', 'vokc'};
    idLayerStereo = {'vopnResponse', 'vokcResponse'};
    nLayer = length(idLayerCorr); % number of layers
    % initialise storage variables
    seedwisePearsonCorrelationPn = cell(nLayer, 1);
    seedwiseDistanceStereotypy = cell(nLayer, 1);
    % find all odor pairs
    nOdor = size(S{1}.simulationData.(idLayerCorr{1}), 2);
    pairOdor = nchoosek(1:nOdor, 2);
    nPairOdor = size(pairOdor, 1);
    % calculate stereotypy values for each layer
    for iLayer = 1:nLayer
        seedwisePearsonCorrelationPn{iLayer} = zeros(nSeed, nPairOdor);
        seedwiseDistanceStereotypy{iLayer} = zeros(nSeed, nPairOdor);
        for iSeed = 1:nSeed
            for iPairOdor = 1:nPairOdor
                seedwisePearsonCorrelationPn{iLayer}(iSeed, iPairOdor) = corr(S{iSeed}.simulationData.(idLayerCorr{iLayer}){1, pairOdor(iPairOdor, 1)}, S{iSeed}.simulationData.(idLayerCorr{iLayer}){1, pairOdor(iPairOdor, 2)});
            end
            seedwiseDistanceStereotypy{iLayer}(iSeed, :) = computeindividualpairstereotypy(spikeData.(idLayerStereo{iLayer}){iSeed});
        end % iSeed
    end % iLayer
    % compute mean, sd and p-value for stereotypy
    % plot the graph for each layer
    objPlot = gramm('y', seedwiseDistanceStereotypy{2}(:), 'x', seedwisePearsonCorrelationPn{1}(:));
    objPlot.geom_point();
    % set plot titles and properties
    objPlot.axe_property('XLim', [-1, 1], 'XGrid', 'on', 'XTick', -1:0.5:1, 'YLim', [-1 1], 'YGrid', 'on', 'YTick', -1:0.5:1);
    objPlot.set_names('y', 'Stereotypy (Total KC response)', 'x', 'Correlation (PN odor vectors)');
    objPlot.set_text_options('font', 'Droid Sans', 'base_size', 8, 'title_scaling', 0.8);
    objPlot.set_title('Total KC response stereotypy vs PN vector correlation for all pairwise odors over all seeds');
    % draw plot, set axis properties and save figure
    handleFigure = figure('Position', [100 100 400 300]);
    objPlot.draw();
    % add correlation values
    text(-0.5, -0.8, computecorrelationstats(seedwiseDistanceStereotypy{2}(:), seedwisePearsonCorrelationPn{1}(:)), 'Parent', objPlot.facet_axes_handles, 'FontSize', 10);
    objPlot.export('file_name', 'plot_stereotypy_vs_correlation', 'export_path', pathPlot, 'file_type', 'png');
    close(handleFigure);
    % save data
    save([pathResult, 'pn_correlation_data.mat'], 'seedwisePearsonCorrelationPn');
    % report successful completion
    run = true;
catch exception
    fprintf('Error saving vector stereotypy vs correlation plots\n')
    txtError = getReport(exception);
    fprintf('%s\n', txtError);
    run = false;
end % try-catch
end % savedistancesquarevectorstereotypyplots
function reportStr = computecorrelationstats(data1, data2)
[r, p] = corr(data1(:), data2(:));
reportStr = sprintf('Pearson''s r: %.4f, P: %1.4g', r, p);
end