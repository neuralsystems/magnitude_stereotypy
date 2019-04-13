function simulateneuronlayer(obj, inputObject)
%%SIMULATENEURONLAYER Simulates the neuron layer and activates all neurons
%whose inputs are greater than the spiking threshold. Modifies the
%rateNeuron property of the neuron layer (Class: VectorSpikingNeuronLayer)
%
% Usage:
%   SIMULATENEURONLAYER(obj, inputObject)
%
% Input:
%    inputObject: object of type NeuronNetwork that encapsulates all
%                 input neuron layers

%**********************************************************************%
% Author: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

% first calculate the total input to each neuron
inputTotal = zeros(obj.nNeuron, 1);
for iInput = 1:obj.nInputLayer
    idInputNeuron = obj.idInputLayer{iInput};
    inputTotal = inputTotal + (inputObject.layerConnectionMatrix.(obj.idInputLayerConnection{iInput}) * inputObject.neuronLayer.(idInputNeuron).rateNeuron);
end
% sets the response as a linear function: slope * input - threshold
obj.rateNeuron = inputTotal * obj.slope - obj.threshold;
end % simulateneuronlayer