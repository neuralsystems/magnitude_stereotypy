function params = defaultparamsvopn()
%%DEFAULTPARAMSVOPN Defines the default parameters for the Vector Olfactory
%Projection Neurons
%
% Usage:
%   params = DEFAULTPARAMSVOPN()
%
% Output:
%   params: structure with all parameters of olfactory projection neurons

%**********************************************************************%
% Programmed and Copyright: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

params = struct('class', @VectorOlfactoryProjectionNeuronLayer, ...
                'idNeuron', 'vopn', ...
                'nNeuron', 50, ...
                'percentNoise', 0, ...
                'rateActiveResponse', [10 30], ...
                'rateBaselineResponse', 0);
end % defaultparamsvopn