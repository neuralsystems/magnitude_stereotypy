classdef VectorOlfactoryProjectionNeuronLayer < handle
    %%VECTOROLFACTORYPROJECTIONNEURONLAYER Inherits from the handle class
    %and implements the odor response of Olfactory Projection Neurons as
    %vectors of firing rates.
    %
    % Usage:
    %   newInstance = VECTOROLFACTORYPROJECTIONNEURONLAYER(params)
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

    properties (SetAccess = immutable)
        idNeuron
        nNeuron
        rateActiveResponse
        rateBaselineResponse
        percentNoise
    end % immutable properties
    %**********************************************************************%
    properties (SetAccess = protected)
        rateNeuron
    end % protected properties
    %**********************************************************************%
    methods
        function obj = VectorOlfactoryProjectionNeuronLayer(params)
            %%VECTOROLFACTORYPROJECTIONNEURONLAYER Class constructor
            obj.idNeuron = params.idNeuron;
            obj.nNeuron = params.nNeuron;
            obj.rateActiveResponse = params.rateActiveResponse;
            obj.rateBaselineResponse = params.rateBaselineResponse;
            obj.percentNoise = params.percentNoise;
        end % class constructor
        %**********************************************************************%
        % declare the implementation of public methods
        resetinputstate(obj, idResponseNeuron, currentSeed, iInput)
    end % methods
    %**********************************************************************%
    methods (Hidden)
        function disp(obj)
            %%DISP Overloads the builtin disp method to display the class
            %information in a specific format
            
            disp(' ');
            stringClassName = sprintf('<a href = "matlab:help %s">%s</a>', class(obj), class(obj));
            disp(['    Input Neuron Layer of Type ', stringClassName]);
            disp(['        Neuron Layer Identifier: ', obj.idNeuron]);
            disp(['        Number of Neurons in the Layer: ',num2str(obj.nNeuron)]);
            stringClassProperties = sprintf('<a href = "matlab:properties %s">properties</a>', class(obj));
            disp(' ');
            disp(['    Display all ', stringClassProperties]);
        end % disp
    end % hidden methods
    %**********************************************************************%
end % VectorOlfactoryProjectionNeuronLayer