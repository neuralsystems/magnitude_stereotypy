classdef VectorOlfactoryProjectionNeuronLayerPartitioned < VectorOlfactoryProjectionNeuronLayer
    %%VECTOROLFACTORYPROJECTIONNEURONLAYERPARTITIONED Inherits from the
    %VectorOlfactoryProjectionNeuronLayer class. The total firing rate for
    %each odor is kept constant and divided across the active neurons
    %
    % Usage:
    %   newInstance = VECTOROLFACTORYPROJECTIONNEURONLAYERPARTITIONED(params)
    %
    % Properties (other than those defined in the parent class):
    %   rateActiveResponseFirst: [min max] range of spiking rate for active neurons for the first odor stimulus
    %
    % See also VECTOROLFACTORYPROJECTIONNEURONLAYER
    
    %**********************************************************************%
    % Author: Aarush Mohit Mittal
    % Contact: aarush (dot) mohit (at) gmail (dot) com
    %**********************************************************************%
    
    properties (SetAccess = immutable)
        % self properties
        rateActiveResponseFirst
    end % immutable properties
    methods
        function obj = VectorOlfactoryProjectionNeuronLayerPartitioned(params)
            %%VECTOROLFACTORYPROJECTIONNEURONLAYERPARTITIONED Class constructor
            obj = obj@VectorOlfactoryProjectionNeuronLayer(params);
            obj.rateActiveResponseFirst = params.rateActiveResponseFirst;
        end % class constructor
        %**********************************************************************%
        % declare the implementation of public methods
        resetinputstate(obj, idResponseNeuron, currentSeed, iInput)
    end % methods
    %**********************************************************************%
end % VectorOlfactoryProjectionNeuronLayerPartitioned