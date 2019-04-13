classdef VectorSpikingNeuronLayer < handle
    %%VECTORSPIKINGNEURONLAYER Implements a Vector Spiking Neuron. It
    %inherits from the handle interface
    %
    % Usage:
    %   newInstance = VECTORSPIKINGNEURONLAYER(params)
    %
    % Inputs:
    %   params: structure of neuron layer properties that are set during initiation of layer instance (should include all immutable properties)
    %
    % Properties (Immutable):
    %                 idNeuron: string specifying the name of the neuron layer
    %                  nNeuron: scalar specifying the number of neurons in the layer
    %             idInputLayer: cell vector of strings specifying the ids of neuron layers that feed into this layer
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
        idNeuron
        nNeuron
        idInputLayer
        nInputLayer
        idInputLayerConnection
        threshold 
    end % immutable properties
    %**********************************************************************%
    properties
        rateNeuron
    end % protected properties
    %**********************************************************************%
    methods
        function obj = VectorSpikingNeuronLayer(params)
            %%VECTORSPIKINGNEURONLAYER Class constructor
            obj.idNeuron = params.idNeuron;
            obj.nNeuron = params.nNeuron;
            obj.idInputLayer = params.idInputLayer;
            obj.nInputLayer = length(obj.idInputLayer);
            inputStringFunction = @(x) [x{:}, '_', obj.idNeuron];
            obj.idInputLayerConnection = arrayfun(inputStringFunction, obj.idInputLayer, 'UniformOutput', false);
            obj.threshold = params.threshold;
        end % class constructor
        %**********************************************************************%
        % declare the implementation of public methods
        simulateneuronlayer(obj, inputObject)
    end % methods
    %**********************************************************************%
    methods (Hidden)
        function disp(obj)
            %%DISP Overloads the builtin disp method to display the class
            %information in a custom format
            
            disp(' ');
            stringClassName = sprintf('<a href = "matlab:help %s">%s</a>', class(obj), class(obj));
            disp(['    Biological Neuron Layer of Type ', stringClassName]);
            disp(['        Neuron Layer Identifier: ', obj.idNeuron]);
            disp(['        Number of Neurons in the Layer: ', num2str(obj.nNeuron)]);
            disp(['        Neuron Layer gets Inputs from: ', strjoin(obj.idInputLayer, ', ')]);
            stringClassProperties = sprintf('<a href = "matlab:properties %s">properties</a>', class(obj));
            disp(' ');
            disp(['    Display all ', stringClassProperties]);
        end % disp
    end % hidden methods
    %**********************************************************************%
end % VectorSpikingNeuronLayeronLayer