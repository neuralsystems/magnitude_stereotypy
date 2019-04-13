function run = savevectorsimulationreports(pathResult, S, nSeed, defaultParams)
%%SAVEVECTORSIMULATIONREPORTS Saves the simulaiton reports for each
%simulation including the parameter report and the spike rate report
%
% Usage:
%   run = SAVEVECTORSIMULATIONREPORTS(pathResult, S, nSeed, defaultParams)
%
% Inputs:
%   defaultParams: struct containing all the default parameters for the simulation
%           nSeed: Number of iterations for each simulation
%               S: matrix of Simulation objects with the simulation data for each iteration
%
% Output:
%   run: logical specifying whether the analysis was run successfully

%**********************************************************************%
% Author: Aarush Mohit Mittal Contact: aarush (dot) mohit (at) gmail (dot)
% com
%**********************************************************************%

try
    pathReport = createresultfolder([pathResult, 'reports\'], true);
    %**********************************************************************%
    % save text spike rate report (averaged over all seeds)
    idNeuronLayer = S{1}.population(1).idNeuronLayer;
    for iSeed = 1:nSeed
        tempRateSpiking = S{iSeed}.rateSpiking.mean;
        for iLayer = idNeuronLayer
            rateSpiking.(iLayer{:}).num(iSeed) = mean(tempRateSpiking.(iLayer{:}).num);
            rateSpiking.(iLayer{:}).all(iSeed) = mean(tempRateSpiking.(iLayer{:}).all);
            rateSpiking.(iLayer{:}).active(iSeed) = mean(tempRateSpiking.(iLayer{:}).active);
            rateSpiking.(iLayer{:}).numodorwise(iSeed, 1:S{iSeed}.nOdor) = accumarray(repelem((1:S{iSeed}.nOdor).', S{iSeed}.simulationParams.nIndividual, 1), tempRateSpiking.(iLayer{:}).num, [], @mean);
            rateSpiking.(iLayer{:}).activeodorwise(iSeed, 1:S{iSeed}.nOdor) = accumarray(repelem((1:S{iSeed}.nOdor).', S{iSeed}.simulationParams.nIndividual, 1), S{iSeed}.rateSpiking.sum.(iLayer{:}).active, [], @mean);
            rateSpiking.(iLayer{:}).baseline(iSeed) = mean(tempRateSpiking.(iLayer{:}).baseline);
        end % for iLayer
    end % for iSeed
    strReport = {sprintf('Spiking Rate Report for Different Neuron Layers (Averaged over %d Seeds)\n', nSeed)};
    strReport = [strReport, sprintf('%%%s%%\n', repmat('-', 1, 50))];
    for iLayer = 1:length(idNeuronLayer)
        strReport = [strReport, sprintf('Number of Active %ss: %.2f ± %.2f\n', convertcase(idNeuronLayer{iLayer}, 'upper'), mean(rateSpiking.(idNeuronLayer{iLayer}).num), std(rateSpiking.(idNeuronLayer{iLayer}).num))];
        strReport = [strReport, sprintf('Spiking Rate of All %ss: %.2f ± %.2f\n', convertcase(idNeuronLayer{iLayer}, 'upper'), mean(rateSpiking.(idNeuronLayer{iLayer}).all), std(rateSpiking.(idNeuronLayer{iLayer}).all))];
        strReport = [strReport, sprintf('Spiking Rate of Active %ss: %.2f ± %.2f\n', convertcase(idNeuronLayer{iLayer}, 'upper'), mean(rateSpiking.(idNeuronLayer{iLayer}).active), std(rateSpiking.(idNeuronLayer{iLayer}).active))];
        strReport = [strReport, sprintf('Spiking Rate of Baseline %ss: %.2f ± %.2f\n', convertcase(idNeuronLayer{iLayer}, 'upper'), mean(rateSpiking.(idNeuronLayer{iLayer}).baseline), std(rateSpiking.(idNeuronLayer{iLayer}).baseline))];
        strReport = [strReport, sprintf('%%%s%%\n', repmat('-', 1, 50))];
    end % for iLayer
    createtextreport(strReport, [pathReport, 'report_spike_rate.txt']);
    save([pathReport, 'spiking_data'], 'rateSpiking');
    %**********************************************************************%
    % save simulation parameter report for first seed
    handleFigure = S{1}.plotsimulationparamsreport(defaultParams, nSeed);
    print(handleFigure, '-dpng', '-r300', [pathReport, sprintf('report_simulation_parameters.png')]);
    close(handleFigure);
    %**********************************************************************%
    run = true;
catch exception
    fprintf('Error saving simulation reports\n')
    errorText = getReport(exception);
    fprintf('%s\n', errorText);
    run = false;
end % try-catch
end % savesimulationreports