function params = defaultparamsvectornetwork()
%%DEFAULTPARAMSVECTORNETWORK Defines the default parameters for the vector neural network
%
% Usage:
%   params = DEFAULTPARAMSVECTORNETWORK()
%
% Output:
%   params: structure with all parameters of neural network

%**********************************************************************%
% Author: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

params.neuronLayer = struct('vopn', defaultparamsvopn(), 'vokc', defaultparamsvokc(), 'vmbon', defaultparamsvmbon());
params.orderLayerSimulation = {'vopn', 'vokc', 'vmbon'};
params.idInputLayer = {'vopn'};
params.matrixVariation = NaN;
params.layerConnection = struct('vopn_vokc', struct('pConnection', [0.14, 0.14], 'fn', 'random', 'vIndividual', true), ...
                                'vokc_vmbon', struct('pConnection', [0.5, 0.5], 'fn', 'non-random-column', 'vIndividual', false));
end % defaultparamsvectornetwork