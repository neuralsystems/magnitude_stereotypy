function params = generatedefaultparams(fnNetworkParams, fnSimulationParams)
%%GENERATEDEFAULTPARAMS Generates the default parameters for both the
%neural network and simulations using the provided functions
%
% Usage:
%   params = GENERATEDEFAULTPARAMS(fnNetworkParams, fnSimulationParams)
%
% Inputs:
%      fnNetworkParams: function defining the network parameters
%   fnSimulationParams: function defining the simulation parameters
%
% Output:
%   params: structure with all parameters

%**********************************************************************%
% Author: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

params.networkParams = fnNetworkParams();
params.simulationParams = fnSimulationParams();
end % generatedefaultparams