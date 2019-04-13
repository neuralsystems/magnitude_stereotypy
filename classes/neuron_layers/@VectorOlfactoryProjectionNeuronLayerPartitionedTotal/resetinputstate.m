function resetinputstate(obj, idResponseNeuron, currentSeed, iInput)
%%RESETINPUTSTATE Sets the state of the neuron layer for the current input
%conditions (Class: VectorOlfactoryProjectionNeuronLayerPartitionedTotal)
%
% Usage:
%   RESETINPUTSTATE(obj, idResponseNeuron, currentSeed, iInput)
%
% Inputs:
%   idResponseNeuron: logical vector specifying the neurons that are active for the current input
%        currentSeed: scalar specifying the random seed for the current stimulus
%             iInput: structure specifying the identity of the current input

%**********************************************************************%
% Author: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

% random number generator is set according to the current input conditions
% for reproducible output
obj.rateNeuron = obj.rateBaselineResponse * ones(obj.nNeuron, 1);
if iInput.(obj.idNeuron) == 1
    sumRateNeuron = round(mean(obj.rateActiveResponse) * obj.nNeuron / 2);
else
    sumRateNeuron = obj.nTotalSpikes;
end
currentSeed = currentSeed + 1e7 * iInput.(obj.idNeuron);
rng(currentSeed, 'twister')
tempRateNeuron = generaterandompartitionmatrixcolumn(sumRateNeuron, sum(idResponseNeuron), obj.rateActiveResponse);
obj.rateNeuron(idResponseNeuron) = tempRateNeuron(randperm(sum(idResponseNeuron)));
end % resetinputstate