function outputMatrix = generaterandommatrix(sizeMatrix, pElement, randomType)
%%GENERATERANDOMMATRIX Generates a logical matrix of specified size using
%the specified method
%
% Usage:
%   outputMatrix = GENERATERANDOMMATRIX(sizeMatrix, pElement, randomType)
%   outputMatrix = GENERATERANDOMMATRIX(sizeMatrix, pElement)
%
% Inputs:
%   sizeMatrix: size of output matrix [nRow nColumn] (1-D and 2-D only)
%     pElement: probability of 'true' elements (number of 'true' elements in case of 'exact-random' methods)
%   randomType: method for randomness generation (Default: 'random')
%               Options:
%                   'random': approximately prod(sizeMatrix) * pElement 'true' elements at random positions
%        'pseudo-random-row': exactly ceil(nColumn * pElement) 'true' elements in each row at random positions
%     'pseudo-random-column': exactly ceil(nRow * pElement) 'true' elements in each column at random positions
%           'non-random-row': first ceil(nRow * pElement) rows will be all 'true'
%        'non-random-column': first ceil(nColumn * pElement) columns will be all 'true'
%         'exact-random-row': exactly pElement 'true' elements in each row at random positions
%      'exact-random-column': exactly pElement 'true' elements in each column at random positions
%
% Output:
%   outputMatrix: logical matrix of specified size
%
% See also RAND, CEIL, PROD, LOGICAL

%**********************************************************************%
% Author: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

% assign default values, and check inputs
if nargin < 3
    randomType = 'random';
end
if length(sizeMatrix) ~= 2
    error('Expected size of output matrix to be 2-D')
end
if ~isscalar(pElement) || pElement < 0
    error('pElement should be a non negative scalar')
end
switch randomType
    %**********************************************************************%
    case 'random'
        if pElement > 1
            error('Probability cannot be greater than 1')
        end
        outputMatrix = logical(rand(sizeMatrix) > (1 - pElement));
    %**********************************************************************%
    case 'pseudo-random-row'
        if pElement > 1
            error('Probability cannot be greater than 1')
        end
        nTrueElement = ceil(pElement * sizeMatrix(2));
        outputMatrix = false(sizeMatrix);
        for iRow = 1:sizeMatrix(1)
            outputMatrix(iRow, randperm(sizeMatrix(2), nTrueElement)) = true;
        end
    %**********************************************************************%    
    case 'pseudo-random-column'
        if pElement > 1
            error('Probability cannot be greater than 1')
        end
        nTrueElement = ceil(pElement * sizeMatrix(1));
        outputMatrix = false(sizeMatrix);
        for iColumn = 1:sizeMatrix(2)
            outputMatrix(randperm(sizeMatrix(1), nTrueElement), iColumn) = true;
        end
    %**********************************************************************%    
    case 'non-random-row'
        if pElement > 1
            error('Probability cannot be greater than 1')
        end
        nTrueElement = ceil(sizeMatrix(1) * pElement);
        outputMatrix = false(sizeMatrix);
        outputMatrix(1:nTrueElement, :) = true;
    %**********************************************************************%
    case 'non-random-column'
        if pElement > 1
            error('Probability cannot be greater than 1')
        end
        nTrueElement = ceil(sizeMatrix(2) * pElement);
        outputMatrix = false(sizeMatrix);
        outputMatrix(:, 1:nTrueElement) = true;
    %**********************************************************************%
    case 'exact-random-row'
        if floor(pElement) ~= pElement
            error('Number of elements should be an integer')
        end
        outputMatrix = false(sizeMatrix);
        for iRow = 1:sizeMatrix(1)
            outputMatrix(iRow, randperm(sizeMatrix(2), pElement)) = true;
        end
    %**********************************************************************%
    case 'exact-random-column'
        if floor(pElement) ~= pElement
            error('Number of elements should be an integer')
        end
        outputMatrix = false(sizeMatrix);
        for iColumn = 1:sizeMatrix(2)
            outputMatrix(randperm(sizeMatrix(1), pElement), iColumn) = true;
        end
    %**********************************************************************%
end
end % generaterandommatrix