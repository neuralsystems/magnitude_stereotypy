classdef VectorOlfactoryProjectionNeuronLayerPartitionedNumber < VectorOlfactoryProjectionNeuronLayer
    %%VECTOROLFACTORYPROJECTIONNEURONLAYERPARTITIONEDNUMBER Inherits from the
    %VectorOlfactoryProjectionNeuronLayer class. The total firing rate for
    %each odor is kept constant and divided across the active neurons. The
    %number of active neurons is different for different odors
    %
    % Usage:
    %   newInstance = VECTOROLFACTORYPROJECTIONNEURONLAYERPARTITIONEDNUMBER(params)
    %
    % Properties (other than those defined in the parent class):
    %   nActiveFirst: the number of active neurons for the first odor stimulus
    %
    % See also VECTOROLFACTORYPROJECTIONNEURONLAYER
    
    %**********************************************************************%
    % Author: Aarush Mohit Mittal
    % Contact: aarush (dot) mohit (at) gmail (dot) com
    %**********************************************************************%
    
    properties (SetAccess = immutable)
        % self properties
        nActiveFirst
    end % immutable properties
    methods
        function obj = VectorOlfactoryProjectionNeuronLayerPartitionedNumber(params)
            %%VECTOROLFACTORYPROJECTIONNEURONLAYERPARTITIONEDNUMBER Class constructor
            obj = obj@VectorOlfactoryProjectionNeuronLayer(params);
            obj.nActiveFirst = params.nActiveFirst;
        end % class constructor
        %**********************************************************************%
        % declare the implementation of public methods
        resetinputstate(obj, idResponseNeuron, currentSeed, iInput)
    end % methods
    %**********************************************************************%
end % VectorOlfactoryProjectionNeuronLayerPartitionedNumber