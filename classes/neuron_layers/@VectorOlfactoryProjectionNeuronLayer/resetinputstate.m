function resetinputstate(obj, idResponseNeuron, currentSeed, iInput)
%%RESETINPUTSTATE Sets the state of the neuron layer for the current input
%conditions (Class: VectorOlfactoryProjectionNeuronLayer)
%
% Usage:
%   RESETINPUTSTATE(obj, idResponseNeuron, currentSeed, iInput)
%
% Inputs:
%   idResponseNeuron: logical vector specifying the neurons that are active for the current stimulus
%        currentSeed: scalar specifying the random seed for the current stimulus
%             iInput: structure specifying the identity of the current stimulus

%**********************************************************************%
% Author: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

% random number generator is set according to the current input conditions
% for reproducible output
obj.rateNeuron = obj.rateBaselineResponse * ones(obj.nNeuron, 1);
newSeed = currentSeed + 1e7 * iInput.(obj.idNeuron);
rng(newSeed, 'twister')
tempRateNeuron = randi(obj.rateActiveResponse, sum(idResponseNeuron), 1);
obj.rateNeuron(idResponseNeuron) = tempRateNeuron(randperm(sum(idResponseNeuron)));
% add gaussian noise across individuals
newSeed = currentSeed + 1e6 * iInput.iIndividual;
rng(newSeed, 'twister')
noise = randn(sum(idResponseNeuron), 1) * obj.percentNoise;
obj.rateNeuron(idResponseNeuron) = max(0, obj.rateNeuron(idResponseNeuron) + noise);
end % resetinputstate