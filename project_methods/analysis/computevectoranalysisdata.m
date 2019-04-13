function computevectoranalysisdata(pathResult, S, nSeed)
%%COMPUTEVECTORANALYSISDATA Converts the raw spike data (for all neuron
%layers) obtained from simulation to analysis ready form and saves it in a
%mat file (spike_data.mat)
%
% Usage:
%   COMPUTEVECTORANALYSISDATA(pathResult, S, nSeed)
%
% Inputs:
%   pathResult: string specifying the current path to the results
%            S: cell matrix containing the simulation data
%        nSeed: scalar specifying the total number of seeds

%**********************************************************************%
% Author: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

% initialise variables for each layer
spikeData.vopnResponse = cell(nSeed, 1);
spikeData.vokcResponse = cell(nSeed, 1);
spikeData.vokcInput = cell(nSeed, 1);
spikeData.vmbonResponse = cell(nSeed, 1);
spikeData.vmbonInput = cell(nSeed, 1);
% calculate responses for each layer
for iSeed = 1:nSeed
    spikeData.vopnResponse{iSeed} = computelayertotalresponse(S{iSeed}, 'vopn');
    spikeData.vokcResponse{iSeed} = computelayertotalresponse(S{iSeed}, 'vokc');
    spikeData.vmbonResponse{iSeed} = computelayertotalresponse(S{iSeed}, 'vmbon');
    tempSpikeData = computetotallayerinput(S{iSeed}, 'vopn', 'vopn_vokc');
    spikeData.vokcInput{iSeed} = cellfun(@sum, tempSpikeData);
    tempSpikeData = computetotallayerinput(S{iSeed}, 'vokc', 'vokc_vmbon');
    spikeData.vmbonInput{iSeed} = cellfun(@sum, tempSpikeData);
end
save([pathResult, 'spike_data.mat'], 'spikeData');
end % computevectoranalysisdata