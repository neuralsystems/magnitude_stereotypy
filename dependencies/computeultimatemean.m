function meanData = computeultimatemean(data)
%%COMPUTEULTIMATEMEAN Finds the mean of the data by averaging over
%all matrix dimensions
%
% Usage:
%   meanData = COMPUTEULTIMATEMEAN(data)
%
% Input:
%   data: n-dimensional matrix
%
% Output:
%   meanData: ultimate mean over all dimensions
%
% Example:
%   data = rand(3,3,3)
%   meanData = COMPUTEULTIMATEMEAN(data) will give a scalar value of mean

%**********************************************************************%
% Author: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

if isempty(data) % ensure mean is not NaN for empty matrices
    meanData = 0;
else % non-empty matrices
    meanData = nanmean(data(:));
end
end % computeultimatemean