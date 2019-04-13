function analysembonstereotypyvectorsimulation(pathResult)
%%ANALYSMBONSTEREOTYPYVECTORSIMULATION Runs analyses for the mbon stereotypy project simulations
%
% Usage:
%   ANALYSMBONSTEREOTYPYVECTORSIMULATION(pathResult)
%
% Input:
%   pathResult: full path to the results folder for the simulation. Should have the generated raw_data.mat file

%**********************************************************************%
% Author: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

if exist([pathResult 'analysis.mat'], 'file')
    load([pathResult 'analysis.mat']);
else
    done = false;
end
if ~done
    % load raw data
    if exist([pathResult, 'raw_data.mat'], 'file')
        load([pathResult, 'raw_data.mat']);
    else
        return
    end
    % calculate analysis data for further use
    if ~exist([pathResult, 'spike_data.mat'], 'file')
        computevectoranalysisdata(pathResult, S, nSeed);
    end
    % save simulation reports
    done = savevectorsimulationreports(pathResult, S, nSeed, defaultParams);
    % calculate stereotypy for various neuron layer inputs and responses
    % save plots for overall statistics and creates text report for mean values
    done = done & saveoverallstereotypyplots(pathResult, nSeed);
    save([pathResult 'analysis.mat'], 'done');
end
end % analysembonstereotypyvectorsimulation