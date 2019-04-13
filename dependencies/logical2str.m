function outputString = logical2str(inputMatrix)
%%LOGICAL2STR Convert logical to string (0: false, 1: true)
%
% Usage:
%   outputString = LOGICAL2STR(inputMatrix)
%
% Inputs:
%   inputMatrix: logical matrix of n-dimension
%
% Output:
%   outputString: String or cell of strings
%
% Example:
%   LOGICAL2STR([true, false])
%   Outputs {'true', 'false}
%
% See also INT2STR, NUM2STR, MAT2STR, FUNC2STR

%**********************************************************************%
% Author: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

if islogical(inputMatrix)
    if isscalar(inputMatrix)
        outputString = strrep([char(inputMatrix * 'true') char(~inputMatrix * 'false')], char(0), '');
    else
        outputString = arrayfun(@logical2str, inputMatrix, 'UniformOutput', false);
    end
else % if isscalar
    error('Expected logical matrix');
end % if islogical
end % logical2str