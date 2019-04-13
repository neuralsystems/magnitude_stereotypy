function run = saveoverallstereotypyplots(pathResult, nSeed)

%**********************************************************************%
% Author: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

try
    % load existing data file
    if exist([pathResult, 'spike_data.mat'], 'file')
        load([pathResult, 'spike_data.mat']);
        if ~exist('spikeData', 'var')
            error('spikeData variable not found in spike_data.mat');
        end
    else
        error('spike_data.mat file doesn''t exist');
    end
    pathPlot = createresultfolder([pathResult, 'stereotypy\'], true);
    % define layer types and titles
    idLayer = {'vopnResponse', 'vokcInput', 'vokcResponse', 'vmbonInput', 'vmbonResponse'};
    titleLayer = {'PN Total Response', 'KC Total Input', 'KC Total Response', 'MBON Total Input', 'MBON Total Response'};
    nLayer = length(idLayer); % number of layers
    % initialise storage variables
    seedwiseDistanceStereotypy = cell(nLayer, 1);
    seedwiseDistanceStereotypyNorm1 = cell(nLayer, 1);
    seedwiseDistanceStereotypyNorm2 = cell(nLayer, 1);
    seedwisePearsonCorrelation = cell(nLayer, 1);
    seedwiseSpearmanCorrelation = cell(nLayer, 1);
    seedwiseCosineSimilarity = cell(nLayer, 1);
    meanStereotypy = cell(nLayer, 1);
    stdStereotypy = cell(nLayer, 1);
    pStereotypy = cell(nLayer, 1);
    % calculate stereotypy values for each layer
    for iLayer = 1:nLayer
        seedwiseDistanceStereotypy{iLayer} = zeros(1, nSeed);
        seedwiseDistanceStereotypyNorm1{iLayer} = zeros(1, nSeed);
        seedwiseDistanceStereotypyNorm2{iLayer} = zeros(1, nSeed);
        seedwisePearsonCorrelation{iLayer} = zeros(1, nSeed);
        seedwiseSpearmanCorrelation{iLayer} = zeros(1, nSeed);
        seedwiseCosineSimilarity{iLayer} = zeros(1, nSeed);
        for iSeed = 1:nSeed
            dataCurrent = spikeData.(idLayer{iLayer}){iSeed};
            seedwiseDistanceStereotypy{iLayer}(iSeed) = computeultimatemean(computeindividualpairstereotypy(dataCurrent));
            seedwiseDistanceStereotypyNorm1{iLayer}(iSeed) = computeultimatemean(computeindividualpairstereotypynormal1(dataCurrent));
            seedwiseDistanceStereotypyNorm2{iLayer}(iSeed) = computeultimatemean(computeindividualpairstereotypynormal2(dataCurrent));
            seedwisePearsonCorrelation{iLayer}(iSeed) = computeultimatemean(computeindividualpairpearsoncorrelation(dataCurrent));
            seedwiseSpearmanCorrelation{iLayer}(iSeed) = computeultimatemean(computeindividualpairspearmancorrelation(dataCurrent));
            seedwiseCosineSimilarity{iLayer}(iSeed) = computeultimatemean(computeindividualpaircosinesimilarity(dataCurrent));
        end % iSeed
        % compute mean, sd and p-value for stereotypy
        meanStereotypy{iLayer} = mean(seedwiseDistanceStereotypy{iLayer});
        stdStereotypy{iLayer} = std(seedwiseDistanceStereotypy{iLayer});
        [~, pStereotypy{iLayer}] = ttest(seedwiseDistanceStereotypy{iLayer}, 0);
        % plot the graph for each layer
        objPlot(1, iLayer) = gramm('x', repmat(titleLayer(iLayer), 1, nSeed), 'y', seedwiseDistanceStereotypy{iLayer}(:));
        objPlot(1, iLayer).geom_jitter('width', 0.5);
        objPlot(1, iLayer).stat_summary('type', 'ci', 'geom', 'black_errorbar');
        % set plot titles and properties
        objPlot(1, iLayer).axe_property('YLim', [-1 1], 'YGrid', 'on', 'YTick', -1:0.5:1);
        objPlot(1, iLayer).set_names('x', '', 'y', 'Stereotypy');
    end % iLayer
    objPlot.set_text_options('font', 'Droid Sans', 'base_size', 8, 'title_scaling', 0.8);
    objPlot.set_title('Overall Stereotypy values for different neuronal layers over all seeds');
    % draw plot, set axis properties and save figure
    handleFigure = figure('Position', [100 100 900 550]);
    objPlot.draw();
    % add mean values
    for iLayer = 1:nLayer
        text(0.5, -0.9, {sprintf('mean: %.4f±%.4f', meanStereotypy{iLayer}, stdStereotypy{iLayer}), sprintf('p-value: %1.4g', pStereotypy{iLayer})}, 'Parent', objPlot(1, iLayer).facet_axes_handles, 'FontSize', 7);
    end % index
    objPlot.export('file_name', 'plot_overall_stereotypy_all_layers', 'export_path', pathPlot, 'file_type', 'png');
    close(handleFigure);
    % save stereotypy values text report
    strReport = {sprintf('Overall Stereotypy Report for Different Neuron Layers (Averaged over %d Seeds)\n\n', nSeed)};
    for iLayer = 1:nLayer
        strReport = [strReport, sprintf('\t%s: %.4f±%.4f, p=%1.4g\n', titleLayer{iLayer}, meanStereotypy{iLayer}, stdStereotypy{iLayer}, pStereotypy{iLayer})];
    end % index
    createtextreport(strReport, [pathResult, 'reports\', 'report_overall_stereotypy.txt']);
    % save data
    save([pathResult, 'overall_stereotypy_data.mat'], 'meanStereotypy', 'stdStereotypy', 'pStereotypy', 'seedwiseDistanceStereotypy', 'seedwiseDistanceStereotypyNorm1', 'seedwiseDistanceStereotypyNorm2', 'seedwisePearsonCorrelation', 'seedwiseSpearmanCorrelation', 'seedwiseCosineSimilarity');
    % report successful completion
    run = true;
catch exception
    fprintf('Error saving vector stereotypy plots and report\n')
    txtError = getReport(exception);
    fprintf('%s\n', txtError);
    run = false;
end % try-catch
end % savedistancesquarevectorstereotypyplots