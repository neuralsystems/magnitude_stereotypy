function rateNeuron = simulateneuronnetwork(obj)
%%SIMULATENEURONNETWORK Simulates the neuron network for the current
%stimulus. All neuron layers are simulated in the order specified
%
% Usage:
%   rateNeuron = SIMULATENEURONNETWORK(obj)
%
% Output:
%   rateNeuron: structure containing the output (spike rates) of all neuron layers

%**********************************************************************%
% Author: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

% resets all neuron layers, generates the spike rate for each neuron layer
% in the order specified
for iNeuron = obj.orderLayerSimulation
    if any(strcmp(iNeuron{:}, obj.idInputLayer))
        continue
    else
        simulateneuronlayer(obj.neuronLayer.(iNeuron{:}), obj);
    end % if idInputLayer
end % for iNeuron
% extract the synaptic output for each neuron layer
for iNeuron = obj.idNeuronLayer
    rateNeuron.(iNeuron{:}) = obj.neuronLayer.(iNeuron{:}).rateNeuron;
end % for iNeuron
end % simulateneuronnetwork