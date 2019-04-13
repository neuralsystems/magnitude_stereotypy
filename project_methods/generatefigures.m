function generatefigures()
%%GENERATEFIGURES Generates all the individual figures for the simulations
%
% Usage:
%   GENERATEFIGURES()

%**********************************************************************%
% Author: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

% set default properties
params.pathCommon = 'mbon_stereotypy_final\';
params.nSeed = 100;
params.pathPlot = createresultfolder([params.pathCommon, 'analyses\']);
params.pathPlotVector = createresultfolder([params.pathCommon, 'analyses\vector\']);
params.valAlpha = 0.5;
params.sizeLine = 2;
params.sizePoint = 5;
params.sizeText = 10;
params.colorMap = [0.2 0.2 0.2];
params.widthErrorBar = 0;
params.widthViolin = 0.3;
params.bandwidthViolin = 0.15;
%---------------------------------------------------------------------------------------------------%

%%%%----bLN1 stereotypy in locust----%%%%
load('data\bln1_stereotypy.mat')
idCategory1 = repelem({'BLN1 response'}, 1, length(corr_values));
idCategory2 = repelem({'BLN1 response'}, 1, length(all_values));
objPlot(1, 1) = gramm('x', idCategory1, 'y', corr_values);
objPlot(1, 1).stat_violin('normalization', 'count', 'fill', 'transparent', 'width', params.widthViolin, 'bandwidth', params.bandwidthViolin);
objPlot(1, 1).set_names('x', '', 'y', 'Correlation stereotypy (locust bLN1 response)');
objPlot(2, 1) = gramm('x', idCategory2, 'y', all_values);
objPlot(2, 1).stat_violin('normalization', 'count', 'fill', 'transparent', 'width', params.widthViolin, 'bandwidth', params.bandwidthViolin);
objPlot(2, 1).set_names('x', '', 'y', 'PRED stereotypy (locust bLN1 response)');
objPlot.axe_property('YLim', [-1 1], 'YGrid', 'on', 'YTick', -1:0.5:1, 'YTickLabel', num2str((-1:0.5:1).', '%.1f'), 'XTickLabel', {}, 'GridLineStyle', '--', 'TickDir', 'out');
objPlot.set_line_options('base_size', params.sizeLine);
objPlot.set_point_options('base_size', params.sizePoint);
objPlot.set_text_options('base_size', params.sizeText);
objPlot.set_color_options('map', [0, 0.663711275482716, 1], 'n_color', 1, 'n_lightness', 1);
% draw plot, set axis properties and save figure
handleFigure = figure('Position', [100, 100, [250 800] + 100]);
rng('default');
objPlot.draw();
% add mean and significance values
text(0.5, -0.8, computeonesamplettest(corr_values, 0), 'Parent', objPlot(1, 1).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(0.5, -0.8, computeonesamplettest(all_values, 0), 'Parent', objPlot(2, 1).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
% export figure
objPlot.export('file_name', 'fig_1bd', 'export_path', params.pathPlot, 'file_type', 'png');
objPlot.export('file_name', 'fig_1bd', 'export_path', params.pathPlotVector, 'file_type', 'pdf');
close(handleFigure);
clearvars -except params
fprintf('Generated Figure 1 (b,d)\n');
%---------------------------------------------------------------------------------------------------%

%%%%----MBON stereotypy in default network----%%%%
pathSimulation = createresultfolder([params.pathCommon, 'default_network_odor_100\']);
load([pathSimulation, 'overall_stereotypy_data.mat']);
idCategory = repelem({'Correlation', 'PRED'}, 1, params.nSeed);
dataPlot = [seedwisePearsonCorrelation{5}, seedwiseDistanceStereotypy{5}];
objPlot(1, 1) = gramm('x', idCategory, 'y', dataPlot);
objPlot(1, 1).stat_violin('normalization', 'count', 'fill', 'transparent', 'width', params.widthViolin * 0.5, 'bandwidth', params.bandwidthViolin);
objPlot(1, 1).set_names('x', ' ', 'y', 'Stereotypy (MBON Response)');
objPlot(1, 1).axe_property('YLim', [0.5 1], 'YGrid', 'on', 'YTick', 0.5:0.1:1, 'YTickLabel', num2str((0.5:0.1:1).', '%.1f'), 'GridLineStyle', '--', 'TickDir', 'out');
objPlot(1, 1).set_color_options('map', [0, 0.663711275482716, 1], 'n_color', 1, 'n_lightness', 1);
objPlot(1, 1).set_layout_options('position', [0 0 0.195 0.5]);
objPlot(1, 1).set_point_options('base_size', params.sizePoint);
%%%%----Pearson correlation vs PRED stereotypy (MBON response) in default network----%%%%
objPlot(1, 2) = gramm('x', seedwiseDistanceStereotypy{5}, 'y', seedwisePearsonCorrelation{5});
objPlot(1, 2).geom_point();
objPlot(1, 2).axe_property('XLim', [0.5 0.9], 'XTick', 0.5:0.1:0.9, 'XTickLabel', num2str((0.5:0.1:0.9).', '%.1f'), 'YLim', [0.95 1], 'YGrid', 'on', 'YTick', 0.95:0.01:1, 'YTickLabel', num2str((0.95:0.01:1).', '%.2f'), 'GridLineStyle', '--', 'TickDir', 'out');
objPlot(1, 2).set_names('x', 'PRED stereotypy (MBON response)', 'y', 'Correlation stereotypy (MBON response)');
objPlot(1, 2).set_color_options('map', [0, 0.663711275482716, 1], 'n_color', 1, 'n_lightness', 1);
objPlot(1, 2).set_layout_options('position', [0.2 0 0.245 0.5]);
objPlot(1, 2).set_point_options('base_size', params.sizePoint);
%%%%----KC stereotypy (individual and overall) in default network----%%%%
total = load([pathSimulation, 'overall_stereotypy_data.mat']);
individual = load([pathSimulation, 'individual_stereotypy_data.mat']);
individual.seedwiseDistanceStereotypy{2}(isnan(individual.seedwiseCorrelation{2})) = [];
individual.seedwiseCorrelation{2}(isnan(individual.seedwiseCorrelation{2})) = [];
% individual
idCategory = repelem({'Correlation', 'PRED'}, 1, length(individual.seedwiseDistanceStereotypy{2}(:)));
dataPlot = [reshape(individual.seedwiseCorrelation{2}, 1, []), reshape(individual.seedwiseDistanceStereotypy{2}, 1, [])];
objPlot(1, 3) = gramm('x', idCategory, 'y', dataPlot);
objPlot(1, 3).stat_violin('normalization', 'count', 'fill', 'transparent', 'width', params.widthViolin, 'bandwidth', params.bandwidthViolin);
objPlot(1, 3).set_names('x', ' ', 'y', 'Stereotypy (Individual KC response)');
objPlot(1, 3).axe_property('YLim', [-1 1], 'YGrid', 'on', 'YTick', -1:0.5:1, 'YTickLabel', num2str((-1:0.5:1).', '%.1f'), 'GridLineStyle', '--', 'TickDir', 'out');
objPlot(1, 3).set_layout_options('position', [0.45 0.5 0.2 0.5]);
objPlot(1, 3).set_point_options('base_size', params.sizePoint);
% overall
idCategory = repelem({'Correlation', 'PRED'}, 1, params.nSeed);
dataPlot = [total.seedwisePearsonCorrelation{3}, total.seedwiseDistanceStereotypy{3}];
objPlot(1, 4) = gramm('x', idCategory, 'y', dataPlot);
objPlot(1, 4).stat_violin('normalization', 'count', 'fill', 'transparent', 'width', params.widthViolin / 2, 'bandwidth', params.bandwidthViolin);
objPlot(1, 4).set_names('x', ' ', 'y', 'Stereotypy (Total KC response)');
objPlot(1, 4).axe_property('YLim', [0.5 1], 'YGrid', 'on', 'YTick', 0.5:0.1:1, 'YTickLabel', num2str((0.5:0.1:1).', '%.1f'), 'GridLineStyle', '--', 'TickDir', 'out');
objPlot(1, 4).set_layout_options('position', [0.45 0 0.2 0.5]);
objPlot(1, 4).set_point_options('base_size', params.sizePoint);
%%%%----across animal correlations in firing rates of all pairs of PNs and KCs----%%%%
if ~exist('data\pn_kc_correlation_data.mat', 'file')
    pathSimulation = createresultfolder([params.pathCommon, 'default_network_odor_100\']);
    load([pathSimulation, 'raw_data.mat']);
    % calculate odor wise distances in firing rates
    nNeuron.pn = S{1}.networkParams.neuronLayer.vopn.nNeuron;
    nNeuron.kc = 50;
    corrCoeff.pn = zeros(nNeuron.pn);
    corrCoeff.kc = zeros(nNeuron.kc);
    [gridX.pn, gridY.pn] = meshgrid(1:nNeuron.pn, 1:nNeuron.pn);
    [gridX.kc, gridY.kc] = meshgrid(1:nNeuron.kc, 1:nNeuron.kc);
    dataRaw = S{1}.rateSpiking.raw;
    pnInd1 = [];
    pnInd2 = [];
    kcInd1 = [];
    kcInd2 = [];
    for iOdor = 1:S{1}.simulationParams.input.vopn.nInput
        pnInd1 = [pnInd1, dataRaw.vopn{1, iOdor}.allNeuron];
        pnInd2 = [pnInd2, dataRaw.vopn{2, iOdor}.allNeuron];
        kcInd1 = [kcInd1, dataRaw.vokc{1, iOdor}.allNeuron];
        kcInd2 = [kcInd2, dataRaw.vokc{2, iOdor}.allNeuron];
    end
    idDel = [];
    for iNeuron = 1:S{1}.networkParams.neuronLayer.vokc.nNeuron
        if all(kcInd1(iNeuron, :) == 0) || all(kcInd2(iNeuron, :) == 0)
            idDel = [idDel, iNeuron];
        end
    end
    kcInd1(idDel, :) = [];
    kcInd2(idDel, :) = [];
    for iNeuron = 1:numel(gridX.pn)
        rho = corr([pnInd1(gridX.pn(iNeuron), :).', pnInd2(gridY.pn(iNeuron), :).']);
        corrCoeff.pn(iNeuron) = rho(2);
    end
    for iNeuron = 1:numel(gridX.kc)
        rho = corr([kcInd1(gridX.kc(iNeuron), :).', kcInd2(gridY.kc(iNeuron), :).']);
        corrCoeff.kc(iNeuron) = rho(2);
    end
    corrCoeff.pn(isnan(corrCoeff.pn)) = 0;
    corrCoeff.kc(isnan(corrCoeff.kc)) = 0;
    save('data\pn_kc_correlation_data.mat', 'gridX', 'gridY', 'corrCoeff')
else
    load('data\pn_kc_correlation_data.mat')
end
% plot correlations
% pn
objPlot(1, 5) = gramm('x', gridX.pn(:), 'y', gridY.pn(:), 'color', corrCoeff.pn(:));
objPlot(1, 5).geom_point();
objPlot(1, 5).axe_property('TickDir', 'out', 'YLim', [1 50], 'YTick', [1 10:10:50], 'XLim', [1 50], 'XTick', [1 10:10:50], 'XTickLabel', num2str([1 10:10:50].', '%.0f'), 'TickDir', 'out');
objPlot(1, 5).set_names('x', 'PNs in individual A', 'y', 'PNs in individual B', 'color', 'correlation');
objPlot(1, 5).set_point_options('markers', {'s'}, 'base_size', params.sizePoint + 1);
objPlot(1, 5).set_continuous_color('CLim', [-1 1]);
objPlot(1, 5).set_layout_options('position', [0.65 0 0.29 0.5]);
objPlot(1, 5).no_legend();
% kc
objPlot(1, 6) = gramm('x', gridX.kc(:), 'y', gridY.kc(:), 'color', corrCoeff.kc(:));
objPlot(1, 6).geom_point();
objPlot(1, 6).axe_property('TickDir', 'out', 'YLim', [1 50], 'YTick', [1 10:10:50], 'XLim', [1 50], 'XTick', [1 10:10:50], 'XTickLabel', num2str([1 10:10:50].', '%.0f'), 'TickDir', 'out');
objPlot(1, 6).set_names('x', 'KCs in individual A', 'y', 'KCs in individual B', 'color', '');
objPlot(1, 6).set_point_options('markers', {'s'}, 'base_size', params.sizePoint + 1);
objPlot(1, 6).set_continuous_color('CLim', [-1 1]);
objPlot(1, 6).set_layout_options('position', [0.65 0.5 0.29 0.5], 'legend_position', [0.95 0.1 0.1 0.9]);
% set common properties, draw plot
objPlot.set_line_options('base_size', params.sizeLine);
objPlot.set_text_options('base_size', params.sizeText);
handleFigure = figure('Position', [100, 100, [1200 600] + 100]);
rng('default');
objPlot.draw();
% add mean and significance values
text(0.5, 0.6, computeonesamplettest(seedwisePearsonCorrelation{5}, 0), 'Parent', objPlot(1, 1).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(1.5, 0.6, computeonesamplettest(seedwiseDistanceStereotypy{5}, 0), 'Parent', objPlot(1, 1).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(0.6, 0.96, computepearsoncorrelation(seedwiseDistanceStereotypy{5}, seedwisePearsonCorrelation{5}), 'Parent', objPlot(1, 2).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(1.5, -0.8, computeonesamplettest(individual.seedwiseDistanceStereotypy{2}, 0), 'Parent', objPlot(1, 3).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(0.5, -0.8, computeonesamplettest(individual.seedwiseCorrelation{2}, 0), 'Parent', objPlot(1, 3).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(1.5, 0.6, computeonesamplettest(total.seedwiseDistanceStereotypy{3}, 0), 'Parent', objPlot(1, 4).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(0.5, 0.6, computeonesamplettest(total.seedwisePearsonCorrelation{3}, 0), 'Parent', objPlot(1, 4).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
% export figure
objPlot.export('file_name', 'fig_2bcdefg', 'export_path', params.pathPlot, 'file_type', 'png');
objPlot.export('file_name', 'fig_2bcdefg', 'export_path', params.pathPlotVector, 'file_type', 'pdf');
close(handleFigure);
clearvars -except params
fprintf('Generated Figure 2 (b,c,d,e,f,g)\n');
%---------------------------------------------------------------------------------------------------%

%%%%----schaffer network cases----%%%%
% load data generated from schaffer 2018 codes
load('data\schaffer_code_data.mat')
%%%%----correlation between PRED stereotypy and correlation stereotypy for untrained network----%%%%
objPlot(1, 1) = gramm('x', dataStereotypyUnmod(:), 'y', dataCorrUnmod(:));
objPlot(1, 1).geom_point('alpha', params.valAlpha);
objPlot(1, 1).set_names('x', 'PRED stereotypy (Piriform readout)', 'y', 'Correlation stereotypy (Piriform readout)');
objPlot(1, 1).axe_property('TickDir', 'out', 'XLim', [-0.02 0.04], 'XTick', -0.02:0.02:0.04, 'XTickLabel', num2str((-0.02:0.02:0.04).', '%.2f'), 'YGrid', 'on', 'YLim', [-0.09 0.08], 'YTick', -0.08:0.04:0.08, 'YTickLabel', num2str((-0.08:0.04:0.08).', '%.2f'), 'GridLineStyle', '--', 'TickDir', 'out');
objPlot(1, 1).set_text_options('base_size', params.sizeText);
objPlot(1, 1).set_color_options('map', params.colorMap, 'n_color', 1, 'n_lightness', 1);
objPlot(1, 1).set_point_options('base_size', params.sizePoint);
objPlot(1, 1).set_layout_options('position', [0 0.5 0.3 0.5], 'legend', false);
%%%%----correlation between PRED stereotypy and correlation stereotypy for trained network----%%%%
objPlot(1, 2) = gramm('x', dataStereotypyUnmodTrained(:), 'y', dataCorrUnmodTrained(:));
objPlot(1, 2).geom_point('alpha', params.valAlpha);
objPlot(1, 2).set_names('x', 'PRED stereotypy (Piriform readout)', 'y', 'Correlation stereotypy (Piriform readout)');
objPlot(1, 2).axe_property('TickDir', 'out', 'XLim', [0.4 0.8], 'XTick', 0.4:0.1:0.8, 'XTickLabel', num2str((0.4:0.1:0.8).', '%.1f'), 'YGrid', 'on', 'YLim', [0.8 1], 'YTick', 0.8:0.04:1, 'YTickLabel', num2str((0.8:0.04:1).', '%.2f'), 'GridLineStyle', '--', 'TickDir', 'out');
objPlot(1, 2).set_text_options('base_size', params.sizeText);
objPlot(1, 2).set_color_options('map', params.colorMap, 'n_color', 1, 'n_lightness', 1);
objPlot(1, 2).set_point_options('base_size', params.sizePoint);
objPlot(1, 2).set_layout_options('position', [0 0 0.3 0.5], 'legend', false);
%%%%----correlation stereotypy with and without weight normalization----%%%%
idColors = repmat({'70%', '30%', '0%'}, 6, 2);
idModification = repelem({'With weight normalization', 'Without weight normalization'}, 6, 3);
dataPlot = [dataCorrUnmod, dataCorrNoMean];
objPlot(1, 3) = gramm('x', idModification(:), 'y', dataPlot(:), 'color', idColors);
objPlot(1, 3).geom_jitter('width', 0.5, 'alpha', params.valAlpha);
objPlot(1, 3).axe_property('YLim', [-0.1, 1], 'YGrid', 'on', 'YTick', 0:0.5:1, 'YTickLabel', num2str((0:0.5:1).', '%.1f'), 'GridLineStyle', '--', 'TickDir', 'out');
objPlot(1, 3).set_names('x', ' ', 'y', 'Correlation stereotypy (Piriform readout)', 'color', '');
objPlot(1, 3).set_order_options('x', 0, 'color', 0);
objPlot(1, 3).set_point_options('base_size', params.sizePoint);
objPlot(1, 3).set_color_options('map', [0.0, 0.7, 0.2; 0, 0.5, 0.6; 0, 0, 0], 'n_color', 3, 'n_lightness', 1);
objPlot(1, 3).set_text_options('base_size', params.sizeText);
objPlot(1, 3).set_layout_options('position', [0.3 0.5 0.3 0.5], 'legend', false);
%%%%----PRED stereotypy with and without weight normalization----%%%%
dataPlot = [dataStereotypyUnmod, dataStereotypyNoMean];
objPlot(1, 4) = gramm('x', idModification(:), 'y', dataPlot(:), 'color', idColors);
objPlot(1, 4).geom_jitter('width', 0.5, 'alpha', params.valAlpha);
objPlot(1, 4).axe_property('YLim', [-0.1, 1], 'YGrid', 'on', 'YTick', 0:0.5:1, 'YTickLabel', num2str((0:0.5:1).', '%.1f'), 'GridLineStyle', '--', 'TickDir', 'out');
objPlot(1, 4).set_names('x', ' ', 'y', 'PRED stereotypy (Piriform readout)', 'color', ' ');
objPlot(1, 4).set_order_options('x', 0, 'color', 0);
objPlot(1, 4).set_point_options('base_size', params.sizePoint);
objPlot(1, 4).set_color_options('map', [0.0, 0.7, 0.2; 0, 0.5, 0.6; 0, 0, 0], 'n_color', 3, 'n_lightness', 1);
objPlot(1, 4).set_text_options('base_size', params.sizeText);
objPlot(1, 4).set_layout_options('position', [0.3 0 0.3 0.5], 'legend_position', [0.55, 0.7, 0.1, 0.1]);
%%%%----default network cases with KC-MBON synapses learning for half the odors randomly----%%%%
learningRange = 0:0.2:1;
stereotypy = zeros(params.nSeed, length(learningRange));
correlation = zeros(params.nSeed, length(learningRange));
idCategory = repmat(reshape(repelem(learningRange, params.nSeed, 1), [], 1), 2, 1);
idMarker = repelem({'PRED'; 'Correlation'}, numel(stereotypy), 1);
for iSim = 1:length(learningRange)
    pathSimulation = createresultfolder([params.pathCommon, sprintf('default_network_odor_100_learn_random_both_p_%.2f\\', learningRange(iSim))]);
    load([pathSimulation, 'overall_stereotypy_data.mat']);
    stereotypy(:, iSim) = seedwiseDistanceStereotypy{5};
    correlation(:, iSim) = seedwisePearsonCorrelation{5};
end
dataPlot = [stereotypy(:); correlation(:)];
objPlot(1, 5) = gramm('x', idCategory, 'y', dataPlot, 'marker', idMarker);
objPlot(1, 5).stat_summary('type', 'sem', 'geom', {'point', 'line', 'black_errorbar'}, 'width', params.widthErrorBar);
objPlot(1, 5).axe_property('YLim', [0 1], 'YGrid', 'on', 'YTick', 0:0.5:1, 'YTickLabel', num2str((0:0.5:1).', '%.1f'), 'XLim', [0 1], 'XTick', 0:0.2:1, 'XTickLabel', num2str((0:0.2:1).', '%.1f'), 'GridLineStyle', '--', 'TickDir', 'out');
objPlot(1, 5).set_names('x', 'Learning rate', 'y', 'Stereotypy (MBON response)');
objPlot(1, 5).set_order_options('marker', 0);
objPlot(1, 5).set_layout_options('position', [0.6 0.5 0.4 0.5], 'legend', false);
objPlot(1, 5).set_line_options('base_size', params.sizeLine);
objPlot(1, 5).set_point_options('base_size', params.sizePoint * 2);
objPlot(1, 5).set_text_options('base_size', params.sizeText);
objPlot(1, 5).set_color_options('map', [0, 0.663711275482716, 1], 'n_color', 1, 'n_lightness', 1);
%%%%----theoretical network with changing KC number----%%%%
if ~exist('data\theoretical_stereotypy.mat', 'file')
    computetheoreticaltotalkcresponsestereotypy;
end
load('data\theoretical_stereotypy.mat');
% plot the stereotypy versus the number of KCs
objPlot(1, 6) = gramm('x', nKRange, 'y', valStereotypy);
objPlot(1, 6).geom_point();
objPlot(1, 6).geom_line();
objPlot(1, 6).axe_property('YLim', [0 1], 'YGrid', 'on', 'YTick', 0:0.5:1, 'YTickLabel', num2str((0:0.5:1).', '%.1f'), 'XLim', [0 1e4], 'XTick', 0:2000:8000, 'XTickLabel', num2str((0:2000:8000).', '%.0f'), 'GridLineStyle', '--', 'TickDir', 'out');
objPlot(1, 6).set_names('x', 'Number of KCs', 'y', 'Expected PRED stereotypy (Total KC response)');
objPlot(1, 6).set_text_options('base_size', params.sizeText);
objPlot(1, 6).set_line_options('base_size', params.sizeLine);
objPlot(1, 6).set_point_options('base_size', params.sizePoint);
objPlot(1, 6).set_layout_options('position', [0.6 0 0.4 0.5]);
% draw plot
handleFigure = figure('Position', [100, 100, [1000 600] + 100]);
rng('default');
objPlot.draw();
% add correlation values
text(-0.01, 0.08, 'Untrained network', 'Parent', objPlot(1, 1).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(0.5, 1, 'Trained network', 'Parent', objPlot(1, 2).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(1, 1, 'Untrained network', 'Parent', objPlot(1, 3).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(1, 1, 'Untrained network', 'Parent', objPlot(1, 4).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(0.01, -0.08, computepearsoncorrelation(dataCorrUnmod(:), dataStereotypyUnmod(:)), 'Parent', objPlot(1, 1).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(0.6, 0.84, computepearsoncorrelation(dataCorrUnmodTrained(:), dataStereotypyUnmodTrained(:)), 'Parent', objPlot(1, 2).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(0.5, 0.2, computeonesamplettest(dataCorrUnmod, 0), 'Parent', objPlot(1, 3).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(1.5, 0.8, computeonesamplettest(dataCorrNoMean, 0), 'Parent', objPlot(1, 3).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(0.5, 0.5, computepairedttest(dataCorrUnmod, dataCorrNoMean), 'Parent', objPlot(1, 3).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(0.5, 0.2, computeonesamplettest(dataStereotypyUnmod, 0), 'Parent', objPlot(1, 4).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(1.5, 0.8, computeonesamplettest(dataStereotypyNoMean, 0), 'Parent', objPlot(1, 4).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(0.5, 0.5, computepairedttest(dataStereotypyUnmod, dataStereotypyNoMean), 'Parent', objPlot(1, 4).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(0.1, 0.1, ['PRED: ', computepearsoncorrelation(learningRange, nanmean(stereotypy, 1))], 'Parent', objPlot(1, 5).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(0.1, 0.15, ['Corr.: ',computepearsoncorrelation(learningRange, nanmean(correlation, 1))], 'Parent', objPlot(1, 5).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
objPlot.export('file_name', 'fig_3', 'export_path', params.pathPlot, 'file_type', 'png');
objPlot.export('file_name', 'fig_3', 'export_path', params.pathPlotVector, 'file_type', 'pdf');
close(handleFigure);
clearvars -except params
fprintf('Generated Figure 3\n');
%---------------------------------------------------------------------------------------------------%

%%%%----total KC input stereotypy in default network----%%%%
pathSimulation = createresultfolder([params.pathCommon, 'default_network\']);
load([pathSimulation, 'overall_stereotypy_data.mat']);
idCategory = repelem({'Total KC input'}, 1, params.nSeed);
objPlot(1, 1) = gramm('x', idCategory, 'y', seedwiseDistanceStereotypy{2});
objPlot(1, 1).stat_violin('normalization', 'count', 'fill', 'transparent', 'width', params.widthViolin / 2, 'bandwidth', params.bandwidthViolin);
objPlot(1, 1).axe_property('YLim', [-1 1], 'YGrid', 'on', 'YTick', -1:0.5:1, 'YTickLabel', num2str((-1:0.5:1).', '%.1f'), 'GridLineStyle', '--', 'TickDir', 'out');
objPlot(1, 1).set_names('x', ' ', 'y', 'Stereotypy (Total KC input)');
objPlot(1, 1).set_color_options('map', [0.028168444512147, 0.714083928404786, 0.294365426432224], 'n_color', 1, 'n_lightness', 1);
objPlot(1, 1).set_layout_options('position', [0 0.5 0.35 0.5]);
%%%%----Generated figure for BLN1 stereotypy in locust within and across concentrations----%%%%
load('data\bln1_stereotypy_within_across.mat');
strCategory = {'Within-concentration', 'Across-concentration'};
idCategory = repelem(strCategory, 1, [length(within_values), length(across_values)]);
plotData = [within_values, across_values];
objPlot(1, 2) = gramm('x', idCategory, 'y', plotData);
objPlot(1, 2).stat_violin('normalization', 'count', 'fill', 'transparent', 'width', 0.5, 'bandwidth', params.bandwidthViolin);
objPlot(1, 2).axe_property('YLim', [-1 1], 'YGrid', 'on', 'YTick', -1:0.5:1, 'YTickLabel', num2str((-1:0.5:1).', '%.1f'), 'GridLineStyle', '--', 'TickDir', 'out');
objPlot(1, 2).set_names('x', ' ', 'y', 'Stereotypy (locust bLN1 response)');
objPlot(1, 2).set_order_options('x', 0);
objPlot(1, 2).set_color_options('map', [0, 0.663711275481716, 1], 'n_color', 1, 'n_lightness', 1);
objPlot(1, 2).set_layout_options('position', [0 0 0.5 0.5]);
%%%%----KC total input/response stereotypy vs difference in PN drive----%%%%
pathSimulation = createresultfolder([params.pathCommon, 'default_network\']);
load([pathSimulation, 'overall_stereotypy_data.mat']);
load([pathSimulation, 'reports\spiking_data.mat']);
rateDiff = abs(diff(rateSpiking.vopn.activeodorwise, 1, 2)).';
%%% KC total input
objPlot(1, 3) = gramm('x', rateDiff, 'y', seedwiseDistanceStereotypy{2}(:));
objPlot(1, 3).geom_point('alpha', params.valAlpha);
objPlot(1, 3).set_names('x', 'Difference in total PN output', 'y', 'Stereotypy (Total KC input)');
objPlot(1, 3).set_color_options('map', [0.028168444512147, 0.714083928404786, 0.294365426432224], 'n_color', 1, 'n_lightness', 1);
objPlot(1, 3).axe_property('YLim', [-1 1], 'YGrid', 'on', 'YTick', -1:0.5:1, 'YTickLabel', num2str((-1:0.5:1).', '%.1f'), 'XLim', [0 330], 'XTick', 0:100:300, 'XTickLabel', num2str((0:100:300).', '%.0f'), 'GridLineStyle', '--', 'TickDir', 'out');
objPlot(1, 3).set_layout_options('position', [0.5 0.5 0.5 0.5]);
%%% KC total response
objPlot(1, 4) = gramm('x', rateDiff, 'y', seedwiseDistanceStereotypy{3}(:));
objPlot(1, 4).geom_point('alpha', params.valAlpha);
objPlot(1, 4).set_names('x', 'Difference in total PN output', 'y', 'Stereotypy (Total KC response)');
objPlot(1, 4).axe_property('YLim', [-1 1], 'YGrid', 'on', 'YTick', -1:0.5:1, 'YTickLabel', num2str((-1:0.5:1).', '%.1f'), 'XLim', [0 330], 'XTick', 0:100:300, 'XTickLabel', num2str((0:100:300).', '%.0f'), 'GridLineStyle', '--', 'TickDir', 'out');
objPlot(1, 4).set_layout_options('position', [0.5 0 0.5 0.5]);
% draw plot, set axis properties and save figure
objPlot.set_line_options('base_size', params.sizeLine);
objPlot.set_point_options('base_size', params.sizePoint);
objPlot.set_text_options('base_size', params.sizeText);
handleFigure = figure('Position', [100, 100, [600 600] + 100]);
rng('default');
objPlot.draw();
% add mean and significance values
text(0.5, -0.8, computeonesamplettest(seedwiseDistanceStereotypy{2}, 0), 'Parent', objPlot(1, 1).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(1, 1, computeunpairedttest(within_values, across_values), 'Parent', objPlot(1, 2).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(0.5, -0.8, computeonesamplettest(within_values, 0), 'Parent', objPlot(1, 2).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(1.5, -0.8, computeonesamplettest(across_values, 0), 'Parent', objPlot(1, 2).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(max(rateDiff) / 3, -0.8, computepearsoncorrelation(rateDiff, seedwiseDistanceStereotypy{2}), 'Parent', objPlot(1, 3).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(max(rateDiff) / 3, -0.8, computepearsoncorrelation(rateDiff, seedwiseDistanceStereotypy{3}), 'Parent', objPlot(1, 4).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
objPlot.export('file_name', 'fig_4', 'export_path', params.pathPlot, 'file_type', 'png');
objPlot.export('file_name', 'fig_4', 'export_path', params.pathPlotVector, 'file_type', 'pdf');
close(handleFigure);
clearvars -except params
fprintf('Generated figure for total KC input stereotypy in default network\n');
%---------------------------------------------------------------------------------------------------%

%%%%----total KC input/response stereotypy in partitioned network----%%%%
% input
pathSimulation = createresultfolder([params.pathCommon, 'partitioned_network\']);
part = load([pathSimulation, 'overall_stereotypy_data.mat']);
objPlot(1, 1) = gramm('x', repelem({'Total KC input'}, 1, params.nSeed), 'y', part.seedwiseDistanceStereotypy{2});
objPlot(1, 1).stat_violin('normalization', 'count', 'fill', 'transparent', 'width', params.widthViolin, 'bandwidth', params.bandwidthViolin);
objPlot(1, 1).set_names('x', ' ', 'y', 'Stereotypy (Total KC input)');
objPlot(1, 1).set_color_options('map', [0.028168444512147, 0.714083928404786, 0.294365426432224], 'n_color', 1, 'n_lightness', 1);
objPlot(1, 1).axe_property('YLim', [-1 1], 'YGrid', 'on', 'YTick', -1:0.5:1, 'YTickLabel', num2str((-1:0.5:1).', '%.1f'), 'GridLineStyle', '--', 'TickDir', 'out');
objPlot(1, 1).set_line_options('base_size', params.sizeLine);
objPlot(1, 1).set_point_options('base_size', params.sizePoint);
objPlot(1, 1).set_text_options('base_size', params.sizeText);
objPlot(1, 1).set_order_options('x', 0);
objPlot(1, 1).set_layout_options('position', [0 0.5 0.2 0.5]);
% response
objPlot(1, 2) = gramm('x', repelem({'Total KC response'}, 1, params.nSeed), 'y', part.seedwiseDistanceStereotypy{3});
objPlot(1, 2).stat_violin('normalization', 'count', 'fill', 'transparent', 'width', params.widthViolin, 'bandwidth', params.bandwidthViolin);
objPlot(1, 2).set_names('x', ' ', 'y', 'Stereotypy (Total KC response)');
objPlot(1, 2).axe_property('YLim', [-1 1], 'YGrid', 'on', 'YTick', -1:0.5:1, 'YTickLabel', num2str((-1:0.5:1).', '%.1f'), 'GridLineStyle', '--', 'TickDir', 'out');
objPlot(1, 2).set_line_options('base_size', params.sizeLine);
objPlot(1, 2).set_point_options('base_size', params.sizePoint);
objPlot(1, 2).set_text_options('base_size', params.sizeText);
objPlot(1, 2).set_order_options('x', 0);
objPlot(1, 2).set_layout_options('position', [0 0 0.2 0.5]);
%%%%----partitioned network cases with first odor at [10 30] and second odor having different firing range----%%%%
rateMin = 10;
rateMax = 30:5:80;
stereotypy.input = zeros(params.nSeed, length(rateMax));
stereotypy.response = zeros(params.nSeed, length(rateMax));
idCategory = cell(params.nSeed, length(rateMax));
for iSim = 1:length(rateMax)
    pathSimulation = createresultfolder([params.pathCommon, sprintf('partitioned_network_both_odor_different_firing_range_%d_%d\\', rateMin, rateMax(iSim))]);
    load([pathSimulation, 'overall_stereotypy_data.mat']);
    stereotypy.input(:, iSim) = seedwiseDistanceStereotypy{2};
    stereotypy.response(:, iSim) = seedwiseDistanceStereotypy{3};
    idCategory(:, iSim) = repelem({sprintf('%d-%d', rateMin, rateMax(iSim))}, params.nSeed, 1);
end
dataPlot = stereotypy.response(:);
idCategoryPlot = idCategory(:);
facetRow = ones(numel(stereotypy.response(:)), 1);
facetCol = ones(numel(stereotypy.response(:)), 1);
%%%%----partitioned network cases with both odors having same firing range----%%%%
rateMin = 10;
rateMax = 30:5:80;
stereotypy.response = zeros(params.nSeed, length(rateMax));
idCategory = cell(params.nSeed, length(rateMax));
for iSim = 1:length(rateMax)
    pathSimulation = createresultfolder([params.pathCommon, sprintf('partitioned_network_both_odor_same_firing_range_%d_%d\\', rateMin, rateMax(iSim))]);
    load([pathSimulation, 'overall_stereotypy_data.mat']);
    stereotypy.response(:, iSim) = seedwiseDistanceStereotypy{3};
    idCategory(:, iSim) = repelem({sprintf('%d-%d', rateMin, rateMax(iSim))}, params.nSeed, 1);
end
dataPlot = [dataPlot; stereotypy.response(:)];
idCategoryPlot = [idCategoryPlot; idCategory(:)];
facetRow = [facetRow; ones(numel(stereotypy.response(:)), 1) * 2];
facetCol = [facetCol; ones(numel(stereotypy.response(:)), 1)];
%%%%----linear partitioned network cases with first odor at [10 30] and second odor having different firing range----%%%%
rateMin = 10;
rateMax = 30:5:80;
stereotypy.response = zeros(params.nSeed, length(rateMax));
idCategory = cell(params.nSeed, length(rateMax));
for iSim = 1:length(rateMax)
    pathSimulation = createresultfolder([params.pathCommon, sprintf('linear_partitioned_network_both_odor_different_firing_range_%d_%d\\', rateMin, rateMax(iSim))]);
    load([pathSimulation, 'overall_stereotypy_data.mat']);
    stereotypy.response(:, iSim) = seedwiseDistanceStereotypy{3};
    idCategory(:, iSim) = repelem({sprintf('%d-%d', rateMin, rateMax(iSim))}, params.nSeed, 1);
end
dataPlot = [dataPlot; stereotypy.response(:)];
idCategoryPlot = [idCategoryPlot; idCategory(:)];
facetRow = [facetRow; ones(numel(stereotypy.response(:)), 1) * 3];
facetCol = [facetCol; ones(numel(stereotypy.response(:)), 1)];
%%%%----partitioned network cases with both odors having different number of active PNs----%%%%
nActivePn = 25:45;
stereotypy.input = zeros(params.nSeed, length(nActivePn));
stereotypy.response = zeros(params.nSeed, length(nActivePn));
idCategory = cell(params.nSeed, length(nActivePn));
for iSim = 1:length(nActivePn)
    pathSimulation = createresultfolder([params.pathCommon, sprintf('partitioned_network_both_odor_different_number_%d\\', nActivePn(iSim))]);
    load([pathSimulation, 'overall_stereotypy_data.mat']);
    stereotypy.input(:, iSim) = seedwiseDistanceStereotypy{2};
    stereotypy.response(:, iSim) = seedwiseDistanceStereotypy{3};
    idCategory(:, iSim) = repelem({sprintf('%d', nActivePn(iSim))}, params.nSeed, 1);
end
dataPlot = [dataPlot; stereotypy.response(:)];
idCategoryPlot = [idCategoryPlot; idCategory(:)];
facetRow = [facetRow; ones(numel(stereotypy.response(:)), 1)];
facetCol = [facetCol; ones(numel(stereotypy.response(:)), 1) * 2];
%%%%----partitioned network cases with both odors having same number of active PNs----%%%%
nActivePn = 25:45;
stereotypy.response = zeros(params.nSeed, length(nActivePn));
idCategory = cell(params.nSeed, length(nActivePn));
for iSim = 1:length(nActivePn)
    pathSimulation = createresultfolder([params.pathCommon, sprintf('partitioned_network_both_odor_same_number_%d\\', nActivePn(iSim))]);
    load([pathSimulation, 'overall_stereotypy_data.mat']);
    stereotypy.response(:, iSim) = seedwiseDistanceStereotypy{3};
    idCategory(:, iSim) = repelem({sprintf('%d', nActivePn(iSim))}, params.nSeed, 1);
end
dataPlot = [dataPlot; stereotypy.response(:)];
idCategoryPlot = [idCategoryPlot; idCategory(:)];
facetRow = [facetRow; ones(numel(stereotypy.response(:)), 1) * 2];
facetCol = [facetCol; ones(numel(stereotypy.response(:)), 1) * 2];
%%%%----partitioned network cases with both odors having different number of active PNs----%%%%
nActivePn = 25:45;
stereotypy.response = zeros(params.nSeed, length(nActivePn));
idCategory = cell(params.nSeed, length(nActivePn));
for iSim = 1:length(nActivePn)
    pathSimulation = createresultfolder([params.pathCommon, sprintf('linear_partitioned_network_both_odor_different_number_%d\\', nActivePn(iSim))]);
    load([pathSimulation, 'overall_stereotypy_data.mat']);
    stereotypy.response(:, iSim) = seedwiseDistanceStereotypy{3};
    idCategory(:, iSim) = repelem({sprintf('%d', nActivePn(iSim))}, params.nSeed, 1);
end
dataPlot = [dataPlot; stereotypy.response(:)];
idCategoryPlot = [idCategoryPlot; idCategory(:)];
facetRow = [facetRow; ones(numel(stereotypy.response(:)), 1) * 3];
facetCol = [facetCol; ones(numel(stereotypy.response(:)), 1) * 2];
% add data to plot
objPlot(1, 3) = gramm('x', idCategoryPlot, 'y', dataPlot);
objPlot(1, 3).facet_grid(facetRow, facetCol, 'scale', 'free_x', 'space', 'fixed', 'row_labels', false, 'column_labels', false);
objPlot(1, 3).stat_summary('type', 'sem', 'geom', {'point', 'black_errorbar'}, 'width', params.widthErrorBar);
objPlot(1, 3).set_names('x', 'cc', 'y', 'Stereotypy (Total KC response)');
objPlot(1, 3).axe_property('YGrid', 'on', 'GridLineStyle', '--', 'TickDir', 'out');
objPlot(1, 3).set_line_options('base_size', params.sizeLine);
objPlot(1, 3).set_point_options('base_size', params.sizePoint * 1.5);
objPlot(1, 3).set_text_options('base_size', params.sizeText);
objPlot(1, 3).set_layout_options('position', [0.2 0 0.5 1]);
%%%%----total KC response stereotypy in shuffled network----%%%%
pathSimulation = createresultfolder([params.pathCommon, 'shuffled_network\']);
shuffled = load([pathSimulation, 'overall_stereotypy_data.mat']);
idCategory = repelem({'Total KC response'}, 1, params.nSeed);
objPlot(1, 4) = gramm('x', idCategory, 'y', shuffled.seedwiseDistanceStereotypy{3});
objPlot(1, 4).stat_violin('normalization', 'count', 'fill', 'transparent', 'width', params.widthViolin, 'bandwidth', params.bandwidthViolin);
objPlot(1, 4).axe_property('YLim', [-1 1], 'YGrid', 'on', 'YTick', -1:0.5:1, 'YTickLabel', num2str((-1:0.5:1).', '%.1f'), 'GridLineStyle', '--', 'TickDir', 'out');
objPlot(1, 4).set_names('x', ' ', 'y', 'Stereotypy (Total KC response)');
objPlot(1, 4).set_line_options('base_size', params.sizeLine);
objPlot(1, 4).set_point_options('base_size', params.sizePoint);
objPlot(1, 4).set_text_options('base_size', params.sizeText);
objPlot(1, 4).set_layout_options('position', [0.7 0.5 0.2 0.5]);
%%%%----theoretical network with changing KC number----%%%%
if ~exist('data\theoretical_stereotypy.mat', 'file')
    computetheoreticaltotalkcresponsestereotypy();
end
load('data\theoretical_stereotypy.mat');
% plot the stereotypy versus the number of KCs
objPlot(1, 5) = gramm('x', nKRange, 'y', zeros(size(nKRange)));
objPlot(1, 5).geom_point();
objPlot(1, 5).geom_line();
objPlot(1, 5).axe_property('YLim', [-0.5 0.5], 'YGrid', 'on', 'YTick', -0.5:0.5:0.5, 'YTickLabel', num2str((-0.5:0.5:0.5).', '%.1f'), 'XLim', [0 1e4], 'XTick', 0:2000:8000, 'XTickLabel', num2str((0:2000:8000).', '%.0f'), 'GridLineStyle', '--', 'TickDir', 'out');
objPlot(1, 5).set_names('x', 'Number of KCs', 'y', 'Expected Stereotypy (Total KC response)');
objPlot(1, 5).set_line_options('base_size', params.sizeLine);
objPlot(1, 5).set_point_options('base_size', params.sizePoint);
objPlot(1, 5).set_text_options('base_size', params.sizeText);
objPlot(1, 5).set_layout_options('position', [0.7 0 0.3 0.5]);
% draw plot
handleFigure = figure('Position', [100, 100, [1000 800] + 100]);
rng('default');
objPlot.draw();
% modify axes for facet plot
set(objPlot(1, 3).facet_axes_handles(1, 1), 'XTick', 1:length(rateMax), 'YLim', [-0.1 1], 'YTick', 0:0.5:1, 'YTickLabel', num2str((0:0.5:1).', '%.1f'));
set(objPlot(1, 3).facet_axes_handles(2, 1), 'XTick', 1:length(rateMax), 'YLim', [-0.1 1], 'YTick', 0:0.5:1, 'YTickLabel', num2str((0:0.5:1).', '%.1f'));
set(objPlot(1, 3).facet_axes_handles(3, 1), 'XTick', 1:length(rateMax), 'XTickLabel', objPlot(1, 3).facet_axes_handles(3, 1).XTickLabel(1:length(rateMax)), 'XTickLabelRotation', 90, 'YLim', [-0.1 1], 'YTick', 0:0.5:1, 'YTickLabel', num2str((0:0.5:1).', '%.1f'));
xlabel(objPlot(1, 3).facet_axes_handles(3, 1), 'Range of PN spiking rates');
set(objPlot(1, 3).facet_axes_handles(1, 2), 'XTick', (length(rateMax) + 1):5:(length(rateMax)+length(nActivePn)), 'YLim', [-0.1 1], 'YTick', 0:0.5:1);
set(objPlot(1, 3).facet_axes_handles(2, 2), 'XTick', (length(rateMax) + 1):5:(length(rateMax)+length(nActivePn)), 'YLim', [-0.1 1], 'YTick', 0:0.5:1);
set(objPlot(1, 3).facet_axes_handles(3, 2), 'XTick', (length(rateMax) + 1):5:(length(rateMax)+length(nActivePn)), 'YLim', [-0.1 1], 'YTick', 0:0.5:1, 'XTickLabel', objPlot(1, 3).facet_axes_handles(3, 2).XTickLabel((length(rateMax) + 1):5:(length(rateMax)+length(nActivePn))));
xlabel(objPlot(1, 3).facet_axes_handles(3, 2), 'Number of active PNs');
% add mean and significance values
text(0.5, -0.8, computeonesamplettest(part.seedwiseDistanceStereotypy{2}, 0), 'Parent', objPlot(1, 1).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(0.5, -0.8, computeonesamplettest(part.seedwiseDistanceStereotypy{3}, 0), 'Parent', objPlot(1, 2).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(0.5, -0.8, computeonesamplettest(shuffled.seedwiseDistanceStereotypy{3}, 0), 'Parent', objPlot(1, 4).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
objPlot.export('file_name', 'fig_5', 'export_path', params.pathPlot, 'file_type', 'png');
objPlot.export('file_name', 'fig_5', 'export_path', params.pathPlotVector, 'file_type', 'pdf');
close(handleFigure);
clearvars -except params
fprintf('Generated Figure 5\n');
%---------------------------------------------------------------------------------------------------%

%%%%----default network cases with changing PN firing rate mean----%%%%
pnMeanVarRange = (-10:5:50) + 20;
stereotypy = zeros(params.nSeed, length(pnMeanVarRange));
number = zeros(params.nSeed, length(pnMeanVarRange));
rate = zeros(params.nSeed, length(pnMeanVarRange));
idCategory = repelem(pnMeanVarRange, params.nSeed, 1);
for iSim = 1:length(pnMeanVarRange)
    pathSimulation = createresultfolder([params.pathCommon, sprintf('default_network_pn_firing_rate_mean_%d\\', pnMeanVarRange(iSim))]);
    load([pathSimulation, 'overall_stereotypy_data.mat']);
    stereotypy(:, iSim) = seedwiseDistanceStereotypy{3};
    load([pathSimulation, '\reports\spiking_data.mat']);
    number(:, iSim) = rateSpiking.vokc.num;
    rate(:, iSim) = rateSpiking.vokc.active;
end
pn.stereotypy = stereotypy;
pn.number = number;
pn.rate = rate;
pn.id = idCategory;
%%%%----default network cases with changing KC number----%%%%
nKcRange = [100 250 500 1000 2e3:2e3:1.8e4];
stereotypy = zeros(params.nSeed, length(nKcRange));
number = zeros(params.nSeed, length(nKcRange));
rate = zeros(params.nSeed, length(nKcRange));
idCategory = repelem(nKcRange, params.nSeed, 1);
for iSim = 1:length(nKcRange)
    pathSimulation = createresultfolder([params.pathCommon, sprintf('default_network_num_kc_%d\\', nKcRange(iSim))]);
    load([pathSimulation, 'overall_stereotypy_data.mat']);
    stereotypy(:, iSim) = seedwiseDistanceStereotypy{3};
    load([pathSimulation, '\reports\spiking_data.mat']);
    number(:, iSim) = rateSpiking.vokc.num;
    rate(:, iSim) = rateSpiking.vokc.active;
end
kc.stereotypy = stereotypy;
kc.number = number;
kc.rate = rate;
kc.id = idCategory;
dataPlot = [pn.stereotypy(:); pn.number(:); pn.rate(:); kc.stereotypy(:); kc.number(:); kc.rate(:)];
idCategoryPlot = [repmat(pn.id(:), 3, 1); repmat(kc.id(:), 3, 1)];
facetRow = [ones(numel(pn.stereotypy), 1); ones(numel(pn.stereotypy), 1) * 2; ones(numel(pn.stereotypy), 1) * 3; ones(numel(kc.stereotypy), 1); ones(numel(kc.stereotypy), 1) * 2; ones(numel(kc.stereotypy), 1) * 3];
facetCol = [ones(numel(pn.stereotypy) * 3, 1); ones(numel(kc.stereotypy) * 3, 1) * 2];
% add data to plot
objPlot(1, 1) = gramm('x', idCategoryPlot, 'y', dataPlot, 'marker', facetCol);
objPlot(1, 1).facet_grid(facetRow, facetCol, 'scale', 'free', 'row_labels', false, 'column_labels', false);
objPlot(1, 1).stat_summary('type', 'sem', 'geom', {'line', 'point', 'black_errorbar'}, 'width', params.widthErrorBar);
objPlot(1, 1).set_names('row', '', 'column', '');
objPlot(1, 1).set_line_options('base_size', params.sizeLine);
objPlot(1, 1).set_point_options('markers', {'o', 's'}, 'base_size', params.sizePoint * 2);
objPlot(1, 1).set_text_options('base_size', params.sizeText);
objPlot(1, 1).axe_property('YGrid', 'on', 'GridLineStyle', '--', 'TickDir', 'out');
objPlot(1, 1).set_color_options('map', params.colorMap, 'n_color', 1, 'n_lightness', 1);
objPlot(1, 1).set_layout_options('position', [0 0 0.6 1], 'legend', false);
%%% stereotypy vs number of active KCs
dataX = [reshape(repmat(mean(pn.number), params.nSeed, 1), [], 1); reshape(repmat(mean(kc.number), params.nSeed, 1), [], 1)];
dataY1 = [pn.stereotypy(:); kc.stereotypy(:)];
dataY2 = [pn.rate(:); kc.rate(:)];
idMarker = [repelem({'Variation in mean spiking rate of PNs'}, length(pn.number(:)), 1); repelem({'Variation in number of KCs'}, length(kc.number(:)), 1)];
objPlot(1, 2) = gramm('x', dataX, 'y', dataY1, 'marker', idMarker);
objPlot(1, 2).stat_summary('type', 'sem', 'geom', {'line', 'point', 'black_errorbar'}, 'width', params.widthErrorBar);
objPlot(1, 2).axe_property('YLim', [0 1], 'YGrid', 'on', 'YTick', 0:0.5:1, 'YTickLabel', num2str((0:0.5:1).', '%.1f'), 'XLim', [0 2e3], 'XTick', 0:500:2e3, 'GridLineStyle', '--', 'TickDir', 'out');
objPlot(1, 2).set_names('x', ' ', 'y', 'Stereotypy (Total KC response)', 'marker', ' ');
objPlot(1, 2).set_line_options('base_size', params.sizeLine);
objPlot(1, 2).set_point_options('base_size', params.sizePoint * 2);
objPlot(1, 2).set_text_options('base_size', params.sizeText);
objPlot(1, 2).set_color_options('map', params.colorMap, 'n_color', 1, 'n_lightness', 1);
objPlot(1, 2).set_layout_options('position', [0.6 0.5 0.4 0.5], 'legend_position', [0.7, 0.48, 0.2, 0.1]);
%%% rate of active KCs vs number of active KCs
objPlot(1, 3) = gramm('x', dataX, 'y', dataY2, 'marker', idMarker);
objPlot(1, 3).stat_summary('type', 'sem', 'geom', {'line', 'point', 'black_errorbar'}, 'width', params.widthErrorBar);
objPlot(1, 3).axe_property('YLim', [0 150], 'YGrid', 'on', 'YTick', 0:50:160, 'XLim', [0 2e3], 'XTick', 0:500:2e3, 'XTickLabel', num2str((0:500:2e3).', '%.0f'), 'GridLineStyle', '--', 'TickDir', 'out');
objPlot(1, 3).set_names('x', 'Number of active KCs', 'y', 'Rate of active KCs', 'marker', '');
objPlot(1, 3).set_line_options('base_size', params.sizeLine);
objPlot(1, 3).set_point_options('base_size', params.sizePoint * 2);
objPlot(1, 3).set_text_options('base_size', params.sizeText);
objPlot(1, 3).set_color_options('map', params.colorMap, 'n_color', 1, 'n_lightness', 1);
objPlot(1, 3).set_layout_options('position', [0.6 0 0.4 0.5], 'legend', false);
% draw plot, set axis properties and save figure
handleFigure = figure('Position', [100, 100, [1000 800] + 100]);
rng('default');
objPlot.draw();
% set axis limits
set(objPlot(1, 1).facet_axes_handles(1, 1), 'YLim', [0 1], 'YTick', 0:0.5:1, 'YTickLabel', num2str((0:0.5:1).', '%.1f'), 'XLim', [10 70]);
set(objPlot(1, 1).facet_axes_handles(1, 2), 'YLim', [0 1], 'YTick', 0:0.5:1, 'XLim', [0 2e4]);
set(objPlot(1, 1).facet_axes_handles(2, 1), 'YLim', [0 2000], 'YTick', 0:500:2000, 'XLim', [10 70]);
set(objPlot(1, 1).facet_axes_handles(2, 2), 'YLim', [0 2000], 'YTick', 0:500:2000, 'XLim', [0 2e4]);
set(objPlot(1, 1).facet_axes_handles(3, 1), 'YLim', [0 150], 'YTick', 0:50:150, 'XLim', [10 70], 'XTickLabel', num2str((10:10:70).', '%.0f'));
set(objPlot(1, 1).facet_axes_handles(3, 2), 'YLim', [0 150], 'YTick', 0:50:150, 'XLim', [0 2e4], 'XTickLabel', num2str((0:5000:2e4).', '%.0f'));
set(objPlot(1, 2).facet_axes_handles, 'XTickLabel', {});
% set axis labels
ylabel(objPlot(1, 1).facet_axes_handles(1, 1), 'Stereotypy (Total KC response)');
ylabel(objPlot(1, 1).facet_axes_handles(2, 1), 'Number of active KCs');
xlabel(objPlot(1, 1).facet_axes_handles(3, 1), 'Mean spiking rate of PNs');
xlabel(objPlot(1, 1).facet_axes_handles(3, 2), 'Number of KCs');
ylabel(objPlot(1, 1).facet_axes_handles(3, 1), 'Spiking rate of active KCs');
% add correlation values
text(20, 0.1, computepearsoncorrelation(pnMeanVarRange, mean(pn.stereotypy)), 'Parent', objPlot(1, 1).facet_axes_handles(1, 1), 'FontSize', params.sizeText * 0.5);
text(20, 100, computepearsoncorrelation(pnMeanVarRange, mean(pn.number)), 'Parent', objPlot(1, 1).facet_axes_handles(2, 1), 'FontSize', params.sizeText * 0.5);
text(20, 10, computepearsoncorrelation(pnMeanVarRange, mean(pn.rate)), 'Parent', objPlot(1, 1).facet_axes_handles(3, 1), 'FontSize', params.sizeText * 0.5);
text(4000, 0.1, computepearsoncorrelation(nKcRange, mean(kc.stereotypy)), 'Parent', objPlot(1, 1).facet_axes_handles(1, 2), 'FontSize', params.sizeText * 0.5);
text(4000, 200, computepearsoncorrelation(nKcRange, mean(kc.number)), 'Parent', objPlot(1, 1).facet_axes_handles(2, 2), 'FontSize', params.sizeText * 0.5);
text(4000, 10, computepearsoncorrelation(nKcRange, mean(kc.rate)), 'Parent', objPlot(1, 1).facet_axes_handles(3, 2), 'FontSize', params.sizeText * 0.5);
objPlot.export('file_name', 'fig_6', 'export_path', params.pathPlot, 'file_type', 'png');
objPlot.export('file_name', 'fig_6', 'export_path', params.pathPlotVector, 'file_type', 'pdf');
close(handleFigure);
clearvars -except params pn kc
fprintf('Generated Figure 6\n');
%---------------------------------------------------------------------------------------------------%

%%%%----BLN1 stereotypy in locust as a function of sparseness for 1sd threshold----%%%%
load('data\bln1_stereotypy.mat')
load('data\bln1_stereotypy_thresh.mat')
strCategory = {'No threshold', 'With threshold'};
idCategory = repelem(strCategory, 1, length(all_values));
plotData = [all_values, all_values_highDep{1}];
objPlot(1, 1) = gramm('x', idCategory, 'y', plotData);
objPlot(1, 1).stat_violin('normalization', 'count', 'fill', 'transparent', 'width', params.widthViolin, 'bandwidth', params.bandwidthViolin);
objPlot(1, 1).axe_property('YLim', [-1 1], 'YGrid', 'on', 'YTick', -1:0.5:1, 'YTickLabel', num2str((-1:0.5:1).', '%.1f'), 'GridLineStyle', '--', 'TickDir', 'out');
objPlot(1, 1).set_names('x', ' ', 'y', 'Stereotypy (locust bLN1 response)');
objPlot(1, 1).set_order_options('x', 0);
objPlot(1, 1).set_line_options('base_size', params.sizeLine);
objPlot(1, 1).set_color_options('map', [0, 0.663711275482716, 1], 'n_color', 1, 'n_lightness', 1);
objPlot(1, 1).set_point_options('base_size', params.sizePoint);
objPlot(1, 1).set_text_options('base_size', params.sizeText);
objPlot(1, 1).set_layout_options('position', [0 0.5 0.5 0.5]);
%%%%----BLN1 stereotypy in locust as a function of sparseness----%%%%
valCategory = 0:4;
idCategory = repelem(0:4, 1, length(all_values));
plotData = [all_values, cell2mat(all_values_highDep)];
objPlot(1, 2) = gramm('x', idCategory, 'y', plotData);
objPlot(1, 2).stat_summary('type', 'sem', 'geom', {'point', 'black_errorbar'}, 'width', params.widthErrorBar);
objPlot(1, 2).stat_glm('geom', {'line'}, 'fullrange', true);
objPlot(1, 2).axe_property('YLim', [-0.1 0.4], 'YGrid', 'on', 'YTick', -0.1:0.1:0.4, 'YTickLabel', num2str((-0.1:0.1:0.4).', '%.1f'), 'XLim', [0 4], 'XTick', valCategory, 'XTickLabel', num2str(valCategory.', '%.0f'), 'GridLineStyle', '--', 'TickDir', 'out');
objPlot(1, 2).set_names('x', 'Threshold (s.d. above mean depolarization)', 'y', 'Stereotypy (locust bLN1 response)');
objPlot(1, 2).set_line_options('base_size', params.sizeLine, 'styles', {'--'});
objPlot(1, 2).set_color_options('map', [0 0 0], 'n_color', 1, 'n_lightness', 1);
objPlot(1, 2).set_point_options('base_size', params.sizePoint);
objPlot(1, 2).set_text_options('base_size', params.sizeText);
objPlot(1, 2).set_layout_options('position', [0.5 0.5 0.5 0.5]);
%%%%----total KC stereotypy in fly w and w/o APL 3 odors----%%%%
if ~exist('data\fly_data.mat', 'file')
    % Import the data
    [~, ~, raw] = xlsread('data\fly_mbon_gcamp3_stereotypy.xlsx', 'stereotypies - 3 odors', 'A5:AD634');
    raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x), raw)) = {''};
    % Replace non-numeric cells with NaN
    R = cellfun(@(x) ~isnumeric(x) && ~islogical(x), raw); % Find non-numeric cells
    raw(R) = {NaN}; % Replace non-numeric cells
    % Create output variable
    data = reshape([raw{:}], size(raw));
    % Allocate imported array to variable names
    control.alpha = reshape(data(:, 13:15), [], 1);
    control.alpha(isnan(control.alpha)) = [];
    apl_tnt.alpha = reshape(data(:, 16:18), [], 1);
    apl_tnt.alpha(isnan(apl_tnt.alpha)) = [];
    control.alpha_p = reshape(data(:, 1:3), [], 1);
    control.alpha_p(isnan(control.alpha_p)) = [];
    apl_tnt.alpha_p = reshape(data(:, 4:6), [], 1);
    apl_tnt.alpha_p(isnan(apl_tnt.alpha_p)) = [];
    control.beta = reshape(data(:, 19:21), [], 1);
    control.beta(isnan(control.beta)) = [];
    apl_tnt.beta = reshape(data(:, 22:24), [], 1);
    apl_tnt.beta(isnan(apl_tnt.beta)) = [];
    control.beta_p = reshape(data(:, 7:9), [], 1);
    control.beta_p(isnan(control.beta_p)) = [];
    apl_tnt.beta_p = reshape(data(:, 10:12), [], 1);
    apl_tnt.beta_p(isnan(apl_tnt.beta_p)) = [];
    control.gamma = reshape(data(:, 25:27), [], 1);
    control.gamma(isnan(control.gamma)) = [];
    apl_tnt.gamma = reshape(data(:, 28:30), [], 1);
    apl_tnt.gamma(isnan(apl_tnt.gamma)) = [];
    save('data\fly_data.mat', 'control', 'apl_tnt')
else
    load('data\fly_data.mat')
end
% plot data
plotData = [cell2mat(struct2cell(control)); cell2mat(struct2cell(apl_tnt))];
idCategory = [repelem({'alpha-lobe', 'alpha''-lobe', 'beta-lobe', 'beta''-lobe', 'gamma-lobe'}, structfun(@length, control)).'; repelem({'alpha-lobe', 'alpha''-lobe', 'beta-lobe', 'beta''-lobe', 'gamma-lobe'}, structfun(@length, apl_tnt)).'];
idLightness = repelem({'Control', 'APL>TNT'}, [sum(structfun(@length, control)), sum(structfun(@length, apl_tnt))]);
objPlot(1, 3) = gramm('x', idCategory, 'y', plotData, 'lightness', idLightness);
objPlot(1, 3).stat_violin('normalization', 'area', 'fill', 'transparent', 'width', params.widthViolin * 2, 'bandwidth', params.bandwidthViolin, 'dodge', 0.6);
objPlot(1, 3).axe_property('YLim', [-1 1], 'YGrid', 'on', 'YTick', -1:0.5:1, 'YTickLabel', num2str((-1:0.5:1).', '%.1f'), 'GridLineStyle', '--', 'TickDir', 'out');
objPlot(1, 3).set_names('x', '', 'y', 'Stereotypy (fly KC population)', 'lightness', ' ');
objPlot(1, 3).set_order_options('x', 0, 'lightness', 0);
objPlot(1, 3).set_color_options('lightness_range', [60 30]);
objPlot(1, 3).set_line_options('base_size', params.sizeLine);
objPlot(1, 3).set_point_options('base_size', params.sizePoint);
objPlot(1, 3).set_text_options('base_size', params.sizeText);
objPlot(1, 3).set_layout_options('position', [0 0 1 0.5], 'legend_position', [0.88, 0.05, 0.1, 0.15]);
% draw plot, set axis properties and save figure
objPlot.set_point_options('base_size', params.sizePoint);
objPlot.set_text_options('base_size', params.sizeText);
handleFigure = figure('Position', [100, 100, [600 400] + 100]);
rng('default');
objPlot.draw();
% add mean and significance values
text(0.5, -0.8, computeonesamplettest(all_values, 0), 'Parent', objPlot(1, 1).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(1.5, -0.8, computeonesamplettest(all_values_highDep{1}, 0), 'Parent', objPlot(1, 1).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(1.25, 1, computepairedttest(all_values_highDep{1}, all_values), 'Parent', objPlot(1, 1).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
all_values_highDep = [{all_values}, all_values_highDep];
text(3, 0.2, computepairwisettest(valCategory, all_values_highDep, 'paired'), 'Parent', objPlot(1, 2).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(0.3, 0.3, computepearsoncorrelation(valCategory.', cellfun(@mean, all_values_highDep).'), 'Parent', objPlot(1, 2).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(0.5, 1, computeunpairedttest(control.alpha, apl_tnt.alpha), 'Parent', objPlot(1, 3).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(0.5, -0.8, computeonesamplettest(control.alpha, 0), 'Parent', objPlot(1, 3).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(0.8, -0.3, computeonesamplettest(apl_tnt.alpha, 0), 'Parent', objPlot(1, 3).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(1.5, 1, computeunpairedttest(control.alpha_p, apl_tnt.alpha_p), 'Parent', objPlot(1, 3).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(1.5, -0.8, computeonesamplettest(control.alpha_p, 0), 'Parent', objPlot(1, 3).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(1.8, -0.3, computeonesamplettest(apl_tnt.alpha_p, 0), 'Parent', objPlot(1, 3).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(2.5, 1, computeunpairedttest(control.beta, apl_tnt.beta), 'Parent', objPlot(1, 3).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(2.5, -0.8, computeonesamplettest(control.beta, 0), 'Parent', objPlot(1, 3).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(2.8, -0.3, computeonesamplettest(apl_tnt.beta, 0), 'Parent', objPlot(1, 3).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(3.5, 1, computeunpairedttest(control.beta_p, apl_tnt.beta_p), 'Parent', objPlot(1, 3).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(3.5, -0.8, computeonesamplettest(control.beta_p, 0), 'Parent', objPlot(1, 3).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(3.8, -0.3, computeonesamplettest(apl_tnt.beta_p, 0), 'Parent', objPlot(1, 3).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(4.5, 1, computeunpairedttest(control.gamma, apl_tnt.gamma), 'Parent', objPlot(1, 3).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(4.5, -0.8, computeonesamplettest(control.gamma, 0), 'Parent', objPlot(1, 3).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(4.8, -0.3, computeonesamplettest(apl_tnt.gamma, 0), 'Parent', objPlot(1, 3).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
% save plot
objPlot.export('file_name', 'fig_7bcd', 'export_path', params.pathPlot, 'file_type', 'png');
objPlot.export('file_name', 'fig_7bcd', 'export_path', params.pathPlotVector, 'file_type', 'pdf');
close(handleFigure);
clearvars -except params
fprintf('Generated Figure 7 (b,c,d)\n');
%---------------------------------------------------------------------------------------------------%

%%%%----default network cases with changing KC-MBON connection probability----%%%%
cKcMbonRange1 = [0.01, 0.02:0.02:0.08, 0.1:0.1:1];
stereotypy1 = zeros(params.nSeed, length(cKcMbonRange1));
idCategory = repelem(cKcMbonRange1, params.nSeed, 1);
for iSim = 1:length(cKcMbonRange1)
    pathSimulation = createresultfolder([params.pathCommon, sprintf('default_network_c_kc_mbon_%.2f\\', cKcMbonRange1(iSim))]);
    load([pathSimulation, 'overall_stereotypy_data.mat']);
    stereotypy1(:, iSim) = seedwiseDistanceStereotypy{5};
end
objPlot(1, 1) = gramm('x', idCategory(:), 'y', stereotypy1(:));
objPlot(1, 1).geom_hline('yintercept', mean(seedwiseDistanceStereotypy{3}), 'style', 'k--');
objPlot(1, 1).stat_summary('type', 'sem', 'geom', {'line', 'point', 'black_errorbar'}, 'width', params.widthErrorBar);
objPlot(1, 1).axe_property('YLim', [0 1], 'YGrid', 'on', 'YTick', 0:0.5:1, 'XLim', [0 1], 'YTickLabel', num2str((0:0.5:1).', '%.1f'), 'XTick', 0:0.2:1, 'XTickLabel', num2str((0:0.2:1).', '%.1f'), 'GridLineStyle', '--', 'TickDir', 'out');
objPlot(1, 1).set_names('x', 'KC-MBON connection probability', 'y', 'Stereotypy (MBON response)');
objPlot(1, 1).set_line_options('base_size', params.sizeLine);
objPlot(1, 1).set_point_options('base_size', params.sizePoint);
objPlot(1, 1).set_text_options('base_size', params.sizeText);
objPlot(1, 1).set_color_options('map', params.colorMap, 'n_color', 1, 'n_lightness', 1);
%%%%----default network cases with percent variation in PN-KC connections across individuals and changing KC-MBON connection probability----%%%%
variationRange = 10 .^ linspace(-2, 0, 21);
cKcMbonRange = 10 .^ linspace(-2, 0, 21);
[V, C] = meshgrid(variationRange, cKcMbonRange);
stereotypy = zeros(size(V));
for iSim = 1:numel(V)
    pathSimulation = createresultfolder([params.pathCommon, sprintf('default_network_c_kc_mbon_%.3f_var_pn_kc_%.3f\\', C(iSim), V(iSim))]);
    load([pathSimulation, 'overall_stereotypy_data.mat']);
    stereotypy(iSim) = mean(seedwiseDistanceStereotypy{5}(:));
end
% perform fit
R = (C ./ V);
myfittype = fittype('x^a/(b+x^a)', 'dependent', {'y'}, 'independent', {'x'}, 'coefficients', {'a', 'b'});
[myfit, gof] = fit(R(:), stereotypy(:), myfittype, 'StartPoint', [0 0]);
fitFn = @(a, b, x) x .^ a ./ (b + (x .^ a));
F = fitFn(myfit.a, myfit.b, R(:));
% plot gradient
objPlot(1, 2) = gramm('x', log10(C(:)), 'y', log10(V(:)), 'color', stereotypy(:));
objPlot(1, 2).geom_point();
objPlot(1, 2).axe_property('YLim', [-2 0], 'YTick', -2:0.5:0, 'YTickLabel', round(10 .^ (-2:0.5:0), 2), 'XLim', [-2 0], 'XTick', -2:0.5:0, 'XTickLabel', round(10 .^ (-2:0.5:0), 2), 'GridLineStyle', '--', 'TickDir', 'out');
objPlot(1, 2).set_names('x', 'KC-MBON connection probability', 'y', 'Randomness in PN-KC connections', 'color', '');
objPlot(1, 2).set_point_options('markers', {'s'}, 'base_size', 13);
objPlot(1, 2).set_text_options('base_size', params.sizeText);
objPlot(1, 2).set_continuous_color('CLim', [0 1]);
% plot fit
objPlot(1, 3) = gramm('x', log10(R(:)), 'y', stereotypy(:));
objPlot(1, 3).stat_summary('type', 'sem', 'geom', {'line'});
objPlot(1, 3).geom_point();
objPlot(1, 3).set_names('x', 'Convergence:randomness ratio', 'y', 'Stereotypy (MBON response)');
objPlot(1, 3).set_order_options('x', 0);
objPlot(1, 3).axe_property('YGrid', 'on', 'YLim', [0 1], 'YTick', 0:0.5:1, 'YTickLabel', num2str((0:0.5:1).', '%.1f'), 'GridLineStyle', '--', 'TickDir', 'out');
objPlot(1, 3).set_text_options('base_size', params.sizeText);
objPlot(1, 3).set_point_options('base_size', params.sizePoint);
objPlot(1, 3).set_line_options('base_size', params.sizeLine);
% draw plot, set axis properties and save figure
objPlot(1, 3).set_color_options('map', [0, 0.663711275482716, 1], 'n_color', 1, 'n_lightness', 1);
handleFigure = figure('Position', [100, 100, [900 200] + 100]);
rng('default');
objPlot.draw();
% add correlation values
text(0.3, 0.1, computepearsoncorrelation(cKcMbonRange1, mean(stereotypy1)), 'Parent', objPlot(1, 1).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
% add fitted line
objPlot(1, 3).update('x', log10(R(:)), 'y', F(:));
objPlot(1, 3).stat_summary('type', 'sem', 'geom', {'line'});
objPlot(1, 3).set_color_options('map', params.colorMap, 'n_color', 1, 'n_lightness', 1);
objPlot.draw();
text(1, 0.2, sprintf('$y=\\frac{x^{%.2f}}{%.2f+x^{%.2f}} $', myfit.a, myfit.b, myfit.a), 'Parent', objPlot(1, 3).facet_axes_handles, 'FontSize', params.sizeText, 'Interpreter', 'latex');
text(1, 0.1, sprintf('$R^2=%.2f$', gof.rsquare), 'Parent', objPlot(1, 3).facet_axes_handles, 'FontSize', params.sizeText * 0.5, 'Interpreter', 'latex');
objPlot.export('file_name', 'fig_8', 'export_path', params.pathPlot, 'file_type', 'png');
objPlot.export('file_name', 'fig_8', 'export_path', params.pathPlotVector, 'file_type', 'pdf');
close(handleFigure);
clearvars -except params
fprintf('Generated Figure 8\n');
%---------------------------------------------------------------------------------------------------%

%%%%----MBON stereotypy with identical connections across individuals in default network----%%%%
pathSimulation = createresultfolder([params.pathCommon, 'default_network_odor_100_same_connection\']);
load([pathSimulation, 'overall_stereotypy_data.mat']);
idCategory = repelem({'Default'}, 1, params.nSeed);
dataPlot = seedwiseDistanceStereotypy{5};
%%%%----MBON stereotypy with identical connections across individuals in default network (with noise)----%%%%
pathSimulation = createresultfolder([params.pathCommon, 'default_network_odor_100_same_connection_with_noise\']);
load([pathSimulation, 'overall_stereotypy_data.mat']);
idCategory = [idCategory, repelem({'With noise'}, 1, params.nSeed)];
dataPlot = [dataPlot seedwiseDistanceStereotypy{5}];
%%%%----MBON stereotypy with identical connections across individuals in default network with odors varying across individuals----%%%%
pathSimulation = createresultfolder([params.pathCommon, 'default_network_odor_100_same_connection_vary_odor_across_individual\']);
load([pathSimulation, 'overall_stereotypy_data.mat']);
idCategory = [idCategory, repelem({'With non-stereotypic PNs'}, 1, params.nSeed)];
dataPlot = [dataPlot seedwiseDistanceStereotypy{5}];
% plot data
objPlot(1, 1) = gramm('x', idCategory, 'y', dataPlot);
objPlot(1, 1).stat_violin('normalization', 'width', 'fill', 'transparent', 'width', params.widthViolin * 2, 'bandwidth', params.bandwidthViolin);
objPlot(1, 1).axe_property('YLim', [-1 1], 'YGrid', 'on', 'YTick', -1:0.5:1, 'YTickLabel', num2str((-1:0.5:1).', '%.1f'), 'GridLineStyle', '--', 'TickDir', 'out');
objPlot(1, 1).set_names('x', '', 'y', 'PRED stereotypy (MBON response)');
objPlot(1, 1).set_line_options('base_size', params.sizeLine);
objPlot(1, 1).set_order_options('x', 0);
objPlot(1, 1).no_legend();
objPlot(1, 1).set_point_options('base_size', params.sizePoint);
objPlot(1, 1).set_text_options('base_size', params.sizeText);
objPlot(1, 1).set_color_options('map', [0, 0.663711275482716, 1], 'n_color', 1, 'n_lightness', 1);
%%%%----default network cases with changing noise levels across individuals----%%%%
noiseRange = 4:4:12;
stereotypy = zeros(params.nSeed, length(noiseRange));
idCategory = repelem(noiseRange, params.nSeed, 1);
for iSim = 1:length(noiseRange)
    pathSimulation = createresultfolder([params.pathCommon, sprintf('default_network_odor_100_with_noise_%d\\', noiseRange(iSim))]);
    load([pathSimulation, 'overall_stereotypy_data.mat']);
    stereotypy(:, iSim) = seedwiseDistanceStereotypy{5};
end
objPlot(2, 1) = gramm('x', idCategory(:), 'y', stereotypy(:));
objPlot(2, 1).stat_violin('normalization', 'width', 'fill', 'transparent', 'width', params.widthViolin * 2, 'bandwidth', params.bandwidthViolin, 'dodge', 1);
objPlot(2, 1).axe_property('YLim', [-1 1], 'YGrid', 'on', 'YTick', -1:0.5:1, 'YTickLabel', num2str((-1:0.5:1).', '%.1f'), 'XLim', [2 noiseRange(end) + 2], 'XTick', noiseRange, 'XTickLabel', num2str(noiseRange.', '%.0f'), 'GridLineStyle', '--', 'TickDir', 'out');
objPlot(2, 1).set_names('x', 'Gaussian noise S.D.', 'y', 'PRED stereotypy (MBON response)');
objPlot(2, 1).set_line_options('base_size', params.sizeLine);
objPlot(2, 1).set_point_options('base_size', params.sizePoint);
objPlot(2, 1).set_text_options('base_size', params.sizeText);
objPlot(2, 1).set_color_options('map', [0, 0.663711275482716, 1], 'n_color', 1, 'n_lightness', 1);
% draw plot, set axis properties and save figure
handleFigure = figure('Position', [100, 100, [400 400] + 100]);
rng('default');
objPlot.draw();
% add mean and significance values
text(0.5, -0.8, computeonesamplettest(dataPlot(1:params.nSeed), 0), 'Parent', objPlot(1, 1).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(1.5, -0.8, computeonesamplettest(dataPlot((params.nSeed + 1):(params.nSeed * 2)), 0), 'Parent', objPlot(1, 1).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(2.5, -0.8, computeonesamplettest(dataPlot((params.nSeed * 2 + 1):end), 0), 'Parent', objPlot(1, 1).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
for iNoise = 1:length(noiseRange)
    text((iNoise - 1) * 4 + 3, -0.8, computeonesamplettest(stereotypy((params.nSeed * (iNoise - 1) + 1):(params.nSeed * iNoise)), 0), 'Parent', objPlot(2, 1).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
end
text(5, 1, computeunpairedttest(stereotypy(1:params.nSeed), stereotypy((params.nSeed + 1):(params.nSeed * 2))), 'Parent', objPlot(2, 1).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(9, 1, computeunpairedttest(stereotypy((params.nSeed + 1):(params.nSeed * 2)), stereotypy((params.nSeed * 2 + 1):end)), 'Parent', objPlot(2, 1).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
objPlot.export('file_name', 'fig_sup_1', 'export_path', params.pathPlot, 'file_type', 'png');
objPlot.export('file_name', 'fig_sup_1', 'export_path', params.pathPlotVector, 'file_type', 'pdf');
close(handleFigure);
clearvars -except params
fprintf('Generated Supplementary Figure 1\n');
%---------------------------------------------------------------------------------------------------%

%%%%----MBON and KC stereotypy in default network comparison across number of individuals----%%%%
% 2 individuals
pathSimulation = createresultfolder([params.pathCommon, 'default_network_ind_2\']);
load([pathSimulation, 'overall_stereotypy_data.mat']);
dataPlot2 = seedwiseDistanceStereotypy{5};
dataPlot1 = seedwisePearsonCorrelation{5};
% 3 individuals
pathSimulation = createresultfolder([params.pathCommon, 'default_network_ind_3\']);
load([pathSimulation, 'overall_stereotypy_data.mat']);
dataPlot2 = [dataPlot2, seedwiseDistanceStereotypy{5}];
dataPlot1 = [dataPlot1, seedwisePearsonCorrelation{5}];
% 4 individuals
pathSimulation = createresultfolder([params.pathCommon, 'default_network_ind_4\']);
load([pathSimulation, 'overall_stereotypy_data.mat']);
dataPlot2 = [dataPlot2, seedwiseDistanceStereotypy{5}];
dataPlot1 = [dataPlot1, seedwisePearsonCorrelation{5}];
% make labels
idInd = [2 3 4];
dataPlot = [dataPlot1, dataPlot2];
idCategory = repelem({'Correlation', 'PRED'}, 1, params.nSeed * length(idInd));
idColor = repmat(repelem(idInd, 1, params.nSeed), 1, 2);
objPlot(1, 1) = gramm('x', idCategory(:), 'y', dataPlot(:), 'lightness', idColor(:));
objPlot(1, 1).stat_violin('normalization', 'count', 'fill', 'transparent', 'width', params.widthViolin * 1.5, 'bandwidth', params.bandwidthViolin, 'dodge', 1);
objPlot(1, 1).set_names('x', 'Number of individuals', 'y', 'Stereotypy (MBON response)', 'lightness', ' ');
objPlot(1, 1).axe_property('YLim', [-1 1], 'YGrid', 'on', 'YTick', -1:0.5:1, 'YTickLabel', num2str((-1:0.5:1).', '%.1f'), 'GridLineStyle', '--', 'TickDir', 'out');
objPlot(1, 1).set_layout_options('legend_position', [0.2, 0.75, 0.15, 0.15]);
objPlot(1, 1).set_color_options('map', [0.504594669032175, 0.893207420837209, 0.922300321608021; 0, 0.669211839158661, 0.725559153932073; 0, 0.445048580214995, 0.537696087497632; 0, 0.232438832967482, 0.359192563403585], 'n_color', 1, 'n_lightness', 4);
strReport1.corr = computepairwisettest(idInd, mat2cell(dataPlot2.', params.nSeed * ones(1, length(idInd)), 1), 'unpaired');
strReport1.pred = computepairwisettest(idInd, mat2cell(dataPlot1.', params.nSeed * ones(1, length(idInd)), 1), 'unpaired');
%%%%----Comparison of the two measures of stereotypy while taking random sets of increasing numbers of odors----%%%%
pathResult = createresultfolder([params.pathCommon, 'default_network_odor_100\']);
if exist([pathResult, 'spike_data.mat'], 'file')
    load([pathResult, 'spike_data.mat']);
    if ~exist('spikeData', 'var')
        error('spikeData variable not found in spike_data.mat');
    end
else
    error('spike_data.mat file doesn''t exist');
end
nIter = 1000;
nGroup = 2:4;
rng(1);
% initialise storage variables
stereotypyPRED = zeros(length(nGroup), nIter);
stereotypyCorr = zeros(length(nGroup), nIter);
% calculate stereotypy values for each layer
for iGroup = nGroup
    for iIter = 1:nIter
        idOdor1 = randperm(100, iGroup);
        dataCurrent = spikeData.vmbonResponse{48}(:, idOdor1);
        stereotypyPRED(iGroup - 1, iIter) = computeultimatemean(computeindividualpairstereotypy(dataCurrent));
        stereotypyCorr(iGroup - 1, iIter) = computeultimatemean(computeindividualpairpearsoncorrelation(dataCurrent));
    end % iIter
end % iGroup
dataPlot = [stereotypyCorr(:); stereotypyPRED(:)];
idColor = repmat(repmat(nGroup.', nIter, 1), 2, 1);
idCategory = repelem({'Correlation'; 'PRED'}, length(nGroup) * nIter, 1);
% plot the graph for each layer
objPlot(2, 1) = gramm('x', idCategory, 'y', dataPlot, 'lightness', idColor);
objPlot(2, 1).stat_violin('normalization', 'count', 'fill', 'transparent', 'width', params.widthViolin * 1.5, 'bandwidth', params.bandwidthViolin, 'dodge', 1);
% set plot titles and properties
objPlot(2, 1).axe_property('YLim', [-1 1], 'YGrid', 'on', 'YTick', -1:0.5:1, 'YTickLabel', num2str((-1:0.5:1).', '%.1f'), 'GridLineStyle', '--', 'TickDir', 'out');
objPlot(2, 1).set_names('x', 'Number of odors', 'y', 'Stereotypy (MBON response)', 'lightness', ' ');
objPlot(2, 1).set_order_options('x', 0, 'lightness', 0);
objPlot(2, 1).set_layout_options('legend_position', [0.2, 0.4, 0.15, 0.15]);
objPlot(2, 1).set_color_options('map', [0.504594669032175, 0.893207420837209, 0.922300321608021; 0, 0.669211839158661, 0.725559153932073; 0, 0.445048580214995, 0.537696087497632; 0, 0.232438832967482, 0.359192563403585], 'n_color', 1, 'n_lightness', 4);
strReport2.corr = computepairwisettest(nGroup, mat2cell(reshape(stereotypyCorr.', [], 1), nIter * ones(1, length(nGroup)), 1), 'unpaired');
strReport2.pred = computepairwisettest(nGroup, mat2cell(reshape(stereotypyPRED.', [], 1), nIter * ones(1, length(nGroup)), 1), 'unpaired');
%%%%----Comparison of the two measures of stereotypy while taking 2 random sets of of 2 odors and their combined set----%%%%
% initialise storage variables
stereotypyPRED = zeros(3, nIter);
stereotypyCorr = zeros(3, nIter);
% calculate stereotypy values for each layer
for iIter = 1:nIter
    idOdor1 = randperm(100, 2);
    idOdor2 = randperm(100, 2);
    idOdorBoth = union(idOdor1, idOdor2);
    while length(idOdorBoth) ~= 4
        idOdor1 = randperm(100, 2);
        idOdor2 = randperm(100, 2);
        idOdorBoth = union(idOdor1, idOdor2);
    end
    dataCurrent = spikeData.vmbonResponse{2}(:, idOdor1);
    stereotypyPRED(1, iIter) = computeultimatemean(computeindividualpairstereotypy(dataCurrent));
    stereotypyCorr(1, iIter) = computeultimatemean(computeindividualpairpearsoncorrelation(dataCurrent));
    dataCurrent = spikeData.vmbonResponse{2}(:, idOdor2);
    stereotypyPRED(2, iIter) = computeultimatemean(computeindividualpairstereotypy(dataCurrent));
    stereotypyCorr(2, iIter) = computeultimatemean(computeindividualpairpearsoncorrelation(dataCurrent));
    dataCurrent = spikeData.vmbonResponse{2}(:, idOdorBoth);
    stereotypyPRED(3, iIter) = computeultimatemean(computeindividualpairstereotypy(dataCurrent));
    stereotypyCorr(3, iIter) = computeultimatemean(computeindividualpairpearsoncorrelation(dataCurrent));
end % iIter
dataPlot = [stereotypyCorr(:); stereotypyPRED(:)];
idColor = repmat(repmat({'1st set'; '2nd set'; 'combined'}, nIter, 1), 2, 1);
idCategory = repelem({'Correlation'; 'PRED'}, 3 * nIter, 1);
% plot the graph for each layer
objPlot(3, 1) = gramm('x', idCategory, 'y', dataPlot, 'lightness', idColor);
objPlot(3, 1).stat_violin('normalization', 'count', 'fill', 'transparent', 'width', params.widthViolin * 1.5, 'bandwidth', params.bandwidthViolin, 'dodge', 1);
% set plot titles and properties
objPlot(3, 1).axe_property('YLim', [-1 1], 'YGrid', 'on', 'YTick', -1:0.5:1, 'YTickLabel', num2str((-1:0.5:1).', '%.1f'), 'GridLineStyle', '--', 'TickDir', 'out');
objPlot(3, 1).set_names('x', '', 'y', 'Stereotypy (MBON response)', 'lightness', ' ');
objPlot(3, 1).set_order_options('x', 0, 'color', 0);
objPlot(3, 1).set_layout_options('legend_position', [0.2, 0.05, 0.15, 0.15]);
objPlot(3, 1).set_color_options('map', [0.504594669032175, 0.893207420837209, 0.922300321608021; 0.504594669032175, 0.893207420837209, 0.922300321608021; 0, 0.445048580214995, 0.537696087497632], 'n_color', 1, 'n_lightness', 3);
strReport3.corr = computepairwisettest([1 2 3], mat2cell(reshape(stereotypyCorr.', [], 1), nIter * ones(1, 3), 1), 'unpaired');
strReport3.pred = computepairwisettest([1 2 3], mat2cell(reshape(stereotypyPRED.', [], 1), nIter * ones(1, 3), 1), 'unpaired');
% draw plot, set common properties and save figure
objPlot.set_line_options('base_size', params.sizeLine);
objPlot.set_point_options('base_size', params.sizePoint);
objPlot.set_text_options('base_size', params.sizeText);
handleFigure = figure('Position', [100, 100, [300 600] + 100]);
rng('default');
objPlot.draw();
% add significance values
text(0.5, 0.75, strReport1.corr, 'Parent', objPlot(1, 1).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(1.5, 0.75, strReport1.pred, 'Parent', objPlot(1, 1).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(0.5, 0.75, strReport2.corr, 'Parent', objPlot(2, 1).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(1.5, 0.75, strReport2.pred, 'Parent', objPlot(2, 1).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(0.5, 0.75, strReport3.corr, 'Parent', objPlot(3, 1).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(1.5, 0.75, strReport3.pred, 'Parent', objPlot(3, 1).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
% export figure
objPlot.export('file_name', 'fig_sup_2bcd', 'export_path', params.pathPlot, 'file_type', 'png');
objPlot.export('file_name', 'fig_sup_2bcd', 'export_path', params.pathPlotVector, 'file_type', 'pdf');
close(handleFigure);
clearvars -except params
fprintf('Generated Supplementary Figure 2 (b,c,d)\n');
%---------------------------------------------------------------------------------------------------%

%%%%----theoretical network distances versus the number of KCs----%%%%
if ~exist('data\theoretical_stereotypy.mat', 'file')
    computetheoreticaltotalkcresponsestereotypy;
end
load('data\theoretical_stereotypy.mat');
dataX = [nKRange, nKRange];
dataY = [log10(D1), log10(D2)];
idColor = [repelem({'E[D_1]'}, 1, length(nKRange)), repelem({'E[D_2]'}, 1, length(nKRange))];
objPlot = gramm('x', dataX, 'y', dataY, 'color', idColor);
objPlot.geom_point();
objPlot.geom_line();
% set plot titles and properties
objPlot.axe_property('YLim', [0 5.5], 'YGrid', 'on', 'YTick', 0:5, 'YTickLabel', num2str((10 .^ (0:5)).', '%.0f'), 'XLim', [0 1e4], 'XTick', 0:2000:8000, 'XTickLabel', num2str((0:2000:8000).', '%.0f'), 'GridLineStyle', '--', 'TickDir', 'out');
objPlot.set_names('x', 'Number of KCs (N_K)', 'y', 'Expected distance (Total KC response)', 'color', '');
objPlot.set_text_options('base_size', params.sizeText);
objPlot.set_line_options('base_size', params.sizeLine);
objPlot.set_point_options('base_size', params.sizePoint);
objPlot.set_layout_options('legend_position', [0.8, 0.7, 0.2, 0.2]);
objPlot.set_color_options('map', [0, 0.663711275482716, 1; 1, 0.367323240931323, 0.413223681757493], 'n_color', 2, 'n_lightness', 1);
% draw plot, set axis properties and save figure
handleFigure = figure('Position', [100, 100, [300 200] + 100]);
rng('default');
objPlot.draw();
objPlot.export('file_name', 'fig_sup_3', 'export_path', params.pathPlot, 'file_type', 'png');
objPlot.export('file_name', 'fig_sup_3', 'export_path', params.pathPlotVector, 'file_type', 'pdf');
close(handleFigure);
clearvars -except params
fprintf('Generated Supplementary Figure 3\n');
%---------------------------------------------------------------------------------------------------%

%%%%----Total KC response PRED stereotypy vs PN response vector correlation----%%%%
pathSimulation = createresultfolder([params.pathCommon, 'default_network_new_seed\']);
if ~exist([pathSimulation, 'pn_correlation_data.mat'], 'file')
    savepnvectorcorrelationvsstereotypyplot(pathSimulation, params.nSeed);
else
    load([pathSimulation, 'pn_correlation_data.mat']);
end
default = load([pathSimulation, 'overall_stereotypy_data.mat'], 'seedwiseDistanceStereotypy');
objPlot(1, 1) = gramm('y', default.seedwiseDistanceStereotypy{3}(:), 'x', seedwisePearsonCorrelationPn{1}(:));
objPlot(1, 1).geom_point();
objPlot(1, 1).axe_property('XLim', [-1 1], 'XTick', -1:0.5:1, 'XTickLabel', num2str((-1:0.5:1).', '%.1f'), 'YLim', [-1 1], 'YGrid', 'on', 'YTick', -1:0.5:1, 'YTickLabel', num2str((-1:0.5:1).', '%.1f'), 'GridLineStyle', '--', 'TickDir', 'out');
objPlot(1, 1).set_names('y', 'Stereotypy (Total KC response)', 'x', 'Correlation in PN response vectors for the two odors)');
%%%%----partitioned network cases with first odor at [10 30] and second odor having different firing range----%%%%
rateMin = 10;
rateMax = 30:5:80;
stereotypy = zeros(params.nSeed, length(rateMax));
idCategory = cell(params.nSeed, length(rateMax));
for iSim = 1:length(rateMax)
    pathSimulation = createresultfolder([params.pathCommon, sprintf('partitioned_network_both_odor_different_firing_range_%d_%d\\', rateMin, rateMax(iSim))]);
    load([pathSimulation, 'overall_stereotypy_data.mat']);
    stereotypy(:, iSim) = seedwiseDistanceStereotypy{2};
    idCategory(:, iSim) = repelem({sprintf('%d-%d', rateMin, rateMax(iSim))}, params.nSeed, 1);
end
objPlot(1, 2) = gramm('x', idCategory(:), 'y', stereotypy(:));
objPlot(1, 2).stat_summary('type', 'sem', 'geom', {'point', 'black_errorbar'}, 'width', params.widthErrorBar);
objPlot(1, 2).axe_property('YLim', [-0.5 0.5], 'YGrid', 'on', 'YTick', -0.5:0.5:0.5, 'YTickLabel', num2str((-0.5:0.5:0.5).', '%.1f'), 'GridLineStyle', '--', 'TickDir', 'out', 'XTickLabelRotation', 90);
objPlot(1, 2).set_names('x', 'Range of PN spiking rates', 'y', 'Stereotypy (Total KC input)');
%%%%----partitioned network cases with both odors having different number of active PNs----%%%%
nActivePn = 25:45;
stereotypy = zeros(params.nSeed, length(nActivePn));
idCategory = repelem(nActivePn, params.nSeed, 1);
for iSim = 1:length(nActivePn)
    pathSimulation = createresultfolder([params.pathCommon, sprintf('partitioned_network_both_odor_different_number_%d\\', nActivePn(iSim))]);
    load([pathSimulation, 'overall_stereotypy_data.mat']);
    stereotypy(:, iSim) = seedwiseDistanceStereotypy{2};
end
% total kc input
objPlot(1, 3) = gramm('x', idCategory(:), 'y', stereotypy(:));
objPlot(1, 3).stat_summary('type', 'sem', 'geom', {'point', 'black_errorbar'}, 'width', params.widthErrorBar);
objPlot(1, 3).axe_property('YLim', [-0.5 0.5], 'YGrid', 'on', 'YTick', -0.5:0.5:0.5, 'YTickLabel', num2str((-0.5:0.5:0.5).', '%.1f'), 'XLim', [25 45], 'XTick', 25:5:45, 'XTickLabel', num2str((25:5:45).', '%.0f'), 'GridLineStyle', '--', 'TickDir', 'out');
objPlot(1, 3).set_names('x', 'Number of active PNs', 'y', 'Stereotypy (Total KC input)');
% draw plot, set axis properties and save figure
objPlot.set_text_options('base_size', params.sizeText);
objPlot.set_line_options('base_size', params.sizeLine);
objPlot.set_point_options('base_size', params.sizePoint);
handleFigure = figure('Position', [100, 100, [900 200] + 100]);
rng('default');
objPlot.draw();
% add correlation values
text(-0.5, -0.8, computepearsoncorrelation(default.seedwiseDistanceStereotypy{3}(:), seedwisePearsonCorrelationPn{1}(:)), 'Parent', objPlot(1, 1).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
objPlot.export('file_name', 'fig_sup_4', 'export_path', params.pathPlot, 'file_type', 'png');
objPlot.export('file_name', 'fig_sup_4', 'export_path', params.pathPlotVector, 'file_type', 'pdf');
close(handleFigure);
clearvars -except params
fprintf('Generated Supplementary Figure 4\n');
%---------------------------------------------------------------------------------------------------%

%%%%----partitioned network cases with both odors having different input drives----%%%%
driveRange = 500:5:540;
stereotypy = zeros(params.nSeed, length(driveRange));
idCategory = repelem(driveRange, params.nSeed, 1);
for iSim = 1:length(driveRange)
    pathSimulation = createresultfolder([params.pathCommon, sprintf('partitioned_network_both_odor_different_input_drive_%d\\', driveRange(iSim))]);
    load([pathSimulation, 'overall_stereotypy_data.mat']);
    stereotypy(:, iSim) = seedwiseDistanceStereotypy{3};
end
dataPlot = stereotypy(:);
idCategoryPlot = idCategory(:);
objPlot = gramm('x', idCategoryPlot, 'y', dataPlot);
objPlot.stat_summary('type', 'sem', 'geom', {'point', 'black_errorbar'}, 'width', params.widthErrorBar);
objPlot.set_names('x', 'Total input drive from PNs', 'y', 'Stereotypy (Total KC response)');
objPlot.axe_property('YLim', [0 1], 'YGrid', 'on', 'YTick', 0:0.5:1, 'YTickLabel', num2str((0:0.5:1).', '%.1f'), 'XLim', [500 540], 'XTickLabel', num2str((500:10:540).', '%.0f'), 'GridLineStyle', '--', 'TickDir', 'out');
objPlot.set_line_options('base_size', params.sizeLine);
objPlot.set_point_options('base_size', params.sizePoint * 2);
objPlot.set_text_options('base_size', params.sizeText);
% draw plot, set axis properties and save figure
handleFigure = figure('Position', [100 100 400 300]);
rng('default');
objPlot.draw();
% save file
objPlot.export('file_name', 'fig_sup_5', 'export_path', params.pathPlot, 'file_type', 'png');
objPlot.export('file_name', 'fig_sup_5', 'export_path', params.pathPlotVector, 'file_type', 'pdf');
close(handleFigure);
clearvars -except params
fprintf('Generated Supplementary Figure 5\n');
%---------------------------------------------------------------------------------------------------%

%%%%----default network cases with changing PN -> KC connection probability----%%%%
cPnKcRange = [0.03 0.05 0.06 0.08 0.1:0.05:1];
stereotypy1 = zeros(params.nSeed, length(cPnKcRange));
idCategory1 = repelem(cPnKcRange, params.nSeed, 1);
for iSim = 1:length(cPnKcRange)
    pathSimulation = createresultfolder([params.pathCommon, sprintf('default_network_c_pn_kc_%.2f\\', cPnKcRange(iSim))]);
    load([pathSimulation, 'overall_stereotypy_data.mat']);
    stereotypy1(:, iSim) = seedwiseDistanceStereotypy{3};
end
objPlot(1, 1) = gramm('x', idCategory1(:), 'y', stereotypy1(:));
objPlot(1, 1).stat_summary('type', 'sem', 'geom', {'line', 'point', 'black_errorbar'}, 'width', params.widthErrorBar);
objPlot(1, 1).axe_property('YLim', [0 1], 'YGrid', 'on', 'YTick', 0:0.5:1, 'YTickLabel', num2str((0:0.5:1).', '%.1f'), 'XLim', [0 1], 'XTick', 0:0.2:1, 'XTickLabel', num2str((0:0.2:1).', '%.1f'), 'GridLineStyle', '--', 'TickDir', 'out');
objPlot(1, 1).set_names('x', 'PN-KC connection probability', 'y', 'Stereotypy (Total KC response)');
%%%%----default network cases with changing KC response threshold----%%%%
thresholdRange = 50:10:250;
stereotypy2 = zeros(params.nSeed, length(thresholdRange));
idCategory2 = repelem(thresholdRange, params.nSeed, 1);
for iSim = 1:length(thresholdRange)
    pathSimulation = createresultfolder([params.pathCommon, sprintf('default_network_kc_thresh_%d\\', thresholdRange(iSim))]);
    load([pathSimulation, 'overall_stereotypy_data.mat']);
    stereotypy2(:, iSim) = seedwiseDistanceStereotypy{3};
end
objPlot(1, 2) = gramm('x', idCategory2(:), 'y', stereotypy2(:));
objPlot(1, 2).stat_summary('type', 'sem', 'geom', {'line', 'point', 'black_errorbar'}, 'width', params.widthErrorBar);
objPlot(1, 2).axe_property('YLim', [0 1], 'YGrid', 'on', 'YTick', 0:0.5:1, 'YTickLabel', num2str((0:0.5:1).', '%.1f'), 'XLim', [50 250], 'XTick', 50:50:250, 'XTickLabel', num2str((50:50:250).', '%.0f'), 'GridLineStyle', '--', 'TickDir', 'out');
objPlot(1, 2).set_names('x', 'KC response threshold', 'y', 'Stereotypy (Total KC response)');
% draw plot, set axis properties and save figure
objPlot.set_line_options('base_size', params.sizeLine);
objPlot.set_point_options('base_size', params.sizePoint);
objPlot.set_text_options('base_size', params.sizeText);
objPlot.set_color_options('map', params.colorMap, 'n_color', 1, 'n_lightness', 1);
handleFigure = figure('Position', [100, 100, [600 200] + 100]);
rng('default');
objPlot.draw();
% add correlation values
text(0.5, 0.1, computepearsoncorrelation(cPnKcRange, mean(stereotypy1)), 'Parent', objPlot(1, 1).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
text(100, 0.1, computepearsoncorrelation(thresholdRange, mean(stereotypy2)), 'Parent', objPlot(1, 2).facet_axes_handles, 'FontSize', params.sizeText * 0.5);
objPlot.export('file_name', 'fig_sup_6', 'export_path', params.pathPlot, 'file_type', 'png');
objPlot.export('file_name', 'fig_sup_6', 'export_path', params.pathPlotVector, 'file_type', 'pdf');
close(handleFigure);
clearvars -except params
fprintf('Generated Supplementary Figure 6\n');
%---------------------------------------------------------------------------------------------------%
end

function reportStr = computeonesamplettest(data, expMean, short)
if nargin < 3
    short = false;
end
data = data(:);
[~, p, ~, stats] = ttest(data, expMean);
if short
    reportStr = sprintf('P: %1.4g', p);
else
    reportStr = {sprintf('mean: %.4f±%.4f vs %.4f', nanmean(data), nanstd(data), expMean), sprintf('t(%d): %.4f, P: %1.4g', stats.df, stats.tstat, p)};
end
end

function reportStr = computepairedttest(data1, data2, short)
if nargin < 3
    short = false;
end
data1 = data1(:);
data2 = data2(:);
[~, p, ~, stats] = ttest(data1, data2);
if short
    reportStr = sprintf('P: %1.4g', p);
else
    reportStr = {sprintf('mean: %.4f±%.4f vs mean: %.4f±%.4f', nanmean(data1), nanstd(data1), nanmean(data2), nanstd(data2)), sprintf('t(%d): %.4f, P: %1.4g', stats.df, stats.tstat, p)};
end
end

function reportStr = computeunpairedttest(data1, data2, short)
if nargin < 3
    short = false;
end
data1 = data1(:);
data2 = data2(:);
[~, p, ~, stats] = ttest2(data1, data2);
if short
    reportStr = sprintf('P: %1.4g', p);
else
    reportStr = sprintf('t(%d): %.4f, P: %1.4g', stats.df, stats.tstat, p);
end
end

function reportStr = computepearsoncorrelation(data1, data2)
[r, p] = corr(data1(:), data2(:));
reportStr = sprintf('Pearson''s r: %.4f, P: %1.4g', r, p);
end

function reportStr = computepairwisettest(idCategory, data, type)
pairData = nchoosek(1:length(idCategory), 2);
nPairData = size(pairData, 1);
reportStr = {};
switch type
    case 'unpaired'
        testFn = @computeunpairedttest;
    case 'paired'
        testFn = @computepairedttest;
end
for iPairData = 1:nPairData
    reportStr = [reportStr, {[sprintf('%d<->%d: ', idCategory(pairData(iPairData, 1)), idCategory(pairData(iPairData, 2))), testFn(data{pairData(iPairData, 1)}, data{pairData(iPairData, 2)}, true)]}];
end
for iData = 1:length(idCategory)
    reportStr = [reportStr, {[sprintf('%d: ', idCategory(iData)), computeonesamplettest(data{iData}, 0, true)]}];
end
end