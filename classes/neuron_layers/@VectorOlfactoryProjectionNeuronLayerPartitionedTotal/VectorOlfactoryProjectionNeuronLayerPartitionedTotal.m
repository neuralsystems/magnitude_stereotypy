classdef VectorOlfactoryProjectionNeuronLayerPartitionedTotal < VectorOlfactoryProjectionNeuronLayer
    %%VECTOROLFACTORYPROJECTIONNEURONLAYERPARTITIONEDTOTAL Inherits from the
    %VectorOlfactoryProjectionNeuronLayer class. The total firing rate for
    %each odor is kept constant and divided across the active neurons. The
    %number of active neurons is different for different odors
    %
    % Usage:
    %   newInstance = VECTOROLFACTORYPROJECTIONNEURONLAYERPARTITIONEDTOTAL(params)
    %
    % Properties (other than those defined in the parent class):
    %   nTotalSpikes: [min max] range of spiking rate for active neurons for the first odor stimulus
    %
    % See also VECTOROLFACTORYPROJECTIONNEURONLAYER
    
    %**********************************************************************%
    % Author: Aarush Mohit Mittal
    % Contact: aarush (dot) mohit (at) gmail (dot) com
    %**********************************************************************%
    
    properties (SetAccess = immutable)
        % self properties
        nTotalSpikes
    end % immutable properties
    methods
        function obj = VectorOlfactoryProjectionNeuronLayerPartitionedTotal(params)
            %%VECTOROLFACTORYPROJECTIONNEURONLAYERPARTITIONEDTOTAL Class constructor
            obj = obj@VectorOlfactoryProjectionNeuronLayer(params);
            obj.nTotalSpikes = params.nTotalSpikes;
        end % class constructor
        %**********************************************************************%
        % declare the implementation of public methods
        resetinputstate(obj, idResponseNeuron, currentSeed, iInput)
    end % methods
    %**********************************************************************%
end % VectorOlfactoryProjectionNeuronLayerPartitionedTotal