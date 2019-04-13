function params = defaultparamsvokc()
%%DEFAULTPARAMSVOKC Defines the default parametrs for the Vector Olfactory
%Kenyon Cells
%
% Usage:
%   params = DEFAULTPARAMSVOKC()
%
% Output:
%   params: structure with all parameters of vector olfactory kenyon cells

%**********************************************************************%
% Programmed and Copyright: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

params = struct('class', @VectorSpikingNeuronLayer, ...
                'idNeuron', 'vokc', ...
                'idInputLayer', {{'vopn'}}, ...
                'nNeuron', 2000, ...
                'threshold', 119);
end % defaultparamsvokc