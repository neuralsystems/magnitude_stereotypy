function params = defaultparamsvmbon()
%%DEFAULTPARAMSVMBON Defines the default parametrs for the Vector Mushroom Body Output Neurons
%
% Usage:
%   params = DEFAULTPARAMSVMBON()
%
% Output:
%   params: structure with all parameters of vector MBONs

%**********************************************************************%
% Programmed and Copyright: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

params = struct('class', @VectorSpikingNeuronLayer, ...
                'idNeuron', 'vmbon', ...
                'idInputLayer', {{'vokc'}}, ...
                'nNeuron', 1, ...
                'threshold', 119);
end % defaultparamsvmbon