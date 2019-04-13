function run = saveindividualstereotypyplots(pathResult, S, nSeed)

%**********************************************************************%
% Author: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

try
    pathPlot = createresultfolder([pathResult, 'stereotypy\'], true);
    % define layer types and titles
    idLayer = {'vopn', 'vokc', 'vmbon'};
    titleLayer = {'PN Total Response', 'KC Total Response', 'MBON Total Response'};
    nLayer = length(idLayer); % number of layers
    % initialise storage variables
    seedwiseDistanceStereotypy = cell(nLayer, 1);
    seedwiseCorrelation = cell(nLayer, 1);
    meanStereotypy = cell(nLayer, 1);
    stdStereotypy = cell(nLayer, 1);
    pStereotypy = cell(nLayer, 1);
    % calculate stereotypy values for each layer
    for iLayer = 1:nLayer
        nNeuron = length(S{1}.simulationData.(idLayer{iLayer}){1});
        seedwiseDistanceStereotypy{iLayer} = zeros(nNeuron, nSeed);
        seedwiseCorrelation{iLayer} = zeros(nNeuron, nSeed);
        for iSeed = 1:nSeed
            dataCurrent{1} = cat(2, S{iSeed}.simulationData.(idLayer{iLayer}){1, :});
            dataCurrent{2} = cat(2, S{iSeed}.simulationData.(idLayer{iLayer}){2, :});
            for iNeuron = 1:nNeuron
                % calculate stereotypy values
                seedwiseDistanceStereotypy{iLayer}(iNeuron, iSeed) = computeultimatemean(computeindividualpairstereotypy([dataCurrent{1}(iNeuron, :); dataCurrent{2}(iNeuron, :)]));
                seedwiseCorrelation{iLayer}(iNeuron, iSeed) = computeindividualpairpearsoncorrelation([dataCurrent{1}(iNeuron, :); dataCurrent{2}(iNeuron, :)]);
            end % iNeuron
        end % iSeed
        % compute mean, sd and p-value for stereotypy
        meanStereotypy{iLayer} = mean(seedwiseDistanceStereotypy{iLayer}(:));
        stdStereotypy{iLayer} = std(seedwiseDistanceStereotypy{iLayer}(:));
        [~, pStereotypy{iLayer}] = ttest(seedwiseDistanceStereotypy{iLayer}(:), 0);
        % plot the graph for each layer
        objPlot(1, iLayer) = gramm('x', repmat(titleLayer(iLayer), 1, nSeed * nNeuron), 'y', seedwiseDistanceStereotypy{iLayer}(:));
        objPlot(1, iLayer).geom_jitter('width', 0.5);
        objPlot(1, iLayer).stat_summary('type', 'ci', 'geom', 'black_errorbar');
        % set plot titles and properties
        objPlot(1, iLayer).axe_property('YLim', [-1 1], 'YGrid', 'on', 'YTick', -1:0.1:1);
        objPlot(1, iLayer).set_names('x', '', 'y', 'Stereotypy');
    end % iLayer
    objPlot.set_text_options('font', 'Droid Sans', 'base_size', 10, 'title_scaling', 0.8);
    objPlot.set_title('Individual Neuron Stereotypy values for different neuronal layers over all seeds');
    % draw plot, set axis properties and save figure
    handleFigure = figure('Position', [100 100 900 550]);
    objPlot.draw();
    % add mean values
    for iLayer = 1:nLayer
        text(1, -0.9, {sprintf('mean: %.4f±%.4f', meanStereotypy{iLayer}, stdStereotypy{iLayer}), sprintf('p-value: %1.4g', pStereotypy{iLayer})}, 'Parent', objPlot(1, iLayer).facet_axes_handles, 'FontSize', 8);
    end % index
    objPlot.export('file_name', 'plot_individual_stereotypy_all_layers', 'export_path', pathPlot, 'file_type', 'png');
    close(handleFigure);
    % save stereotypy values text report
    strReport = {sprintf('Individual Neuron Stereotypy Report for Different Neuron Layers (Averaged over %d Seeds)\n\n', nSeed)};
    for iLayer = 1:nLayer
        strReport = [strReport, sprintf('\t%s: %.4f±%.4f, p=%1.4g\n', titleLayer{iLayer}, meanStereotypy{iLayer}, stdStereotypy{iLayer}, pStereotypy{iLayer})];
    end % index
    createtextreport(strReport, [pathResult, 'reports\', 'report_individual_stereotypy.txt']);
    % save data
    save([pathResult, 'individual_stereotypy_data.mat'], 'meanStereotypy', 'stdStereotypy', 'pStereotypy', 'seedwiseDistanceStereotypy', 'seedwiseCorrelation');
    % report successful completion
    run = true;
catch exception
    fprintf('Error saving vector stereotypy plots and report\n')
    txtError = getReport(exception);
    fprintf('%s\n', txtError);
    run = false;
end % try-catch
end % savedistancesquarevectorstereotypyplots