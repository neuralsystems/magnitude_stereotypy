classdef VectorOlfactoryProjectionNeuronLayerShuffled < VectorOlfactoryProjectionNeuronLayer
    %%VECTOROLFACTORYPROJECTIONNEURONLAYERSHUFFLED Inherits from the
    %VectorOlfactoryProjectionNeuronLayer abstract class and implements
    %the odor response for Olfactory Projection Neurons as vectors of
    %firing rates. Each odor is a shuffled version of the other
    %
    % Usage:
    %   newInstance = VectorOlfactoryProjectionNeuronLayerShuffled(params)
    %
    % Inputs:
    %   params: structure of neuron layer properties that are set during initiation of layer instance (should include all immutable properties)
    %
    % Properties (Immutable):
    %               idNeuron: string specifying the name of the neuron layer
    %                nNeuron: scalar specifying the number of neurons in the layer
    %     rateActiveResponse: [min max] range of spiking rate for active neurons for each input
    %   rateBaselineResponse: scalar specifying the baseline spiking rate for neurons
    %           percentNoise: scalar specifying the standard deviation of the gaussian noise across individuals
    %
    % Properties (Calculated):
    %             rateNeuron: spike rate of each neuron (size: nNeuron x 1)
    %
    % Methods:
    %        resetinputstate: resets the state of the neuron layer before simulation of each new input
    
    %**********************************************************************%
    % Author: Aarush Mohit Mittal
    % Contact: aarush (dot) mohit (at) gmail (dot) com
    %**********************************************************************%
    
    methods
        function obj = VectorOlfactoryProjectionNeuronLayerShuffled(params)
            %%VECTOROLFACTORYPROJECTIONNEURONLAYERSHUFFLED Class constructor
            obj = obj@VectorOlfactoryProjectionNeuronLayer(params);
        end % class constructor
        %**********************************************************************%
        % declare the implementation of public methods
        resetinputstate(obj, idResponseNeuron, currentSeed, iInput)
    end % methods
    %**********************************************************************%
end % VectorOlfactoryProjectionNeuronLayerShuffled