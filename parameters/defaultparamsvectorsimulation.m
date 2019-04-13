function params = defaultparamsvectorsimulation()
%%DEFAULTPARAMSVECTORSIMULATION Defines the default parameters for the vector simulations
%
% Usage:
%   params = DEFAULTPARAMSVECTORSIMULATION()
%
% Output:
%   params: structure with all parameters of the simulations

%**********************************************************************%
% Author: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

params.nIndividual = 2;
params.input = struct('vopn', struct('nInput', 2, 'pInput', 0.5, 'vInput', struct('individual', false, 'input', true), 'fnInput', 'random'));
end % defaultparamsvectorsimulation