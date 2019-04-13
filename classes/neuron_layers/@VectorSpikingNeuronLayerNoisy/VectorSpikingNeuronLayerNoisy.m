classdef VectorSpikingNeuronLayerNoisy < VectorSpikingNeuronLayer
    %%VECTORSPIKINGNEURONLAYERNOISY Implements a Vector Spiking Neuron. It
    %inherits from the VectorSpikingNeuronLayer Class and adds noise to all
    %inputs
    %
    % Usage:
    %   newInstance = VECTORSPIKINGNEURONLAYERNOISY(params)
    %
    % Inputs:
    %   params: structure of neuron layer properties that are set during initiation of layer instance (should include all immutable properties)
    %
    % Properties (Immutable):
    %       idNeuron: string specifying the name of the neuron layer
    %        nNeuron: scalar specifying the number of neurons in the layer
    %   idInputLayer: cell vector of strings specifying the ids of neuron layers that feed into this layer
    %   percentNoise: scalar specifying the standard deviation of the gaussian noise across individuals
    %
    % Properties (Calculated):
    %              nInputLayer: scalar specifying the number of input layers
    %   idInputLayerConnection: cell vector of strings specifying the ids of connections from input layers
    %               rateNeuron: vector that stores the spike rate of each neuron (size: nNeuron x 1)
    %                threshold: scalar specifying the spiking threshold for each neuron
    %
    % Methods:
    %    simulateneuronlayer: simulates the neuron layer for the current stimulus
    
    %**********************************************************************%
    % Programmed and Copyright: Aarush Mohit Mittal
    % Contact: aarush (dot) mohit (at) gmail (dot) com
    %**********************************************************************%

    properties (SetAccess = immutable)
        percentNoise
    end % immutable properties
    %**********************************************************************%
    methods
        function obj = VectorSpikingNeuronLayerNoisy(params)
            %%VECTORSPIKINGNEURONLAYERNOISY Class constructor
            obj = obj@VectorSpikingNeuronLayer(params);
            obj.percentNoise = params.percentNoise;
        end % class constructor
        %**********************************************************************%
        % declare the implementation of public methods
        simulateneuronlayer(obj, inputObject)
    end % methods
    %**********************************************************************%
end % VectorSpikingNeuronLayeronLayerNoisy