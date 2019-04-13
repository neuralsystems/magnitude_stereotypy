function outputMatrix = generaterandompartitionmatrixcolumn(sumColumn, nRow, limits, method)
%%GENERATERANDOMPARTITIONMATRIXCOLUMN Generates a random matrix wherein the sum
%of columns is fixed (fixed margins)
%
% Usage:
%   outputMatrix = GENERATERANDOMPARTITIONMATRIXCOLUMN(sumColumn, nRow, limits, method)
%   outputMatrix = GENERATERANDOMPARTITIONMATRIXCOLUMN(sumColumn, nRow, limits)
%   outputMatrix = GENERATERANDOMPARTITIONMATRIXCOLUMN(sumColumn, nRow)
%
% Inputs:
%    sumColumn: row matrix of size 1 x n denoting the sum of each column of the resulting matrix
%         nRow: scalar denoting the number of rows in the resulting matrix
%       limits: [minimum maximum] value allowed in the matrix (Default: [0 max(sumColumn)])
%       method: method for partition calculations (Default: 0)
%
% Output:
%   outputMatrix: matrix of size m x n with fixed margins as defined above

%**********************************************************************%
% Author: Aarush Mohit Mittal
% Contact: aarush (dot) mohit (at) gmail (dot) com
%**********************************************************************%

if ~isrow(sumColumn)
    warning('Sum of columns is not a row vector. Reshaping it...')
    sumColumn = reshape(sumColumn, 1, []);
end
if nargin < 4
    method = 0;
end
if nargin < 3 || length(limits) > 2 || length(limits) < 1
    limits = [0, max(sumColumn)];
end
if length(limits) == 1
    limits = [limits, limits];
end
minVal = limits(1);
maxVal = limits(2);
if any(sumColumn > (maxVal * nRow))
    warning('Maximum value not compatible with given column sums. Setting maximum value to max(sumColumn)...')
    maxVal = max(sumColumn);
end
if any(sumColumn < (minVal * nRow))
    warning('Minimum value not compatible with given column sums. Setting minimum value to 0...')
    minVal = 0;
end
nColumn = length(sumColumn);
outputMatrix = ones(nRow, nColumn) * minVal;
newSumColumn = sumColumn - nRow * minVal;
newMaxVal = maxVal - minVal;
% generates each column by partitioning the desired column sum into nRow
% integers.
if method == 1
    nGroups = floor(newSumColumn / newMaxVal);
    for iColumn = 1:nColumn
        tempColumn = sum(generaterandommatrix([nRow, newMaxVal], nGroups(iColumn), 'exact-random-column'), 2);
        sumLeft = newSumColumn(iColumn) - sum(tempColumn);
        if sumLeft > 0
            nLeft = sum(tempColumn < newMaxVal);
            tempColumn(tempColumn < newMaxVal) = tempColumn(tempColumn < newMaxVal) + generaterandommatrix([nLeft, 1], sumLeft, 'exact-random-column');
        end
        outputMatrix(:, iColumn) = outputMatrix(:, iColumn) + tempColumn;
    end % for iColumn
elseif method == 2
    for iColumn = 1:nColumn
        tempColumn = zeros(nRow, newMaxVal);
        tempColumn(randperm(numel(tempColumn), newSumColumn(iColumn))) = 1;
        tempColumn = sum(tempColumn, 2);
        sumLeft = newSumColumn(iColumn) - sum(tempColumn);
        if sumLeft > 0
            nLeft = sum(tempColumn < newMaxVal);
            tempColumn(tempColumn < newMaxVal) = tempColumn(tempColumn < newMaxVal) + generaterandommatrix([nLeft, 1], sumLeft, 'exact-random-column');
        end
        outputMatrix(:, iColumn) = outputMatrix(:, iColumn) + tempColumn;
    end % for iColumn
elseif method == 3
    for iColumn = 1:nColumn
        tempColumn = randi([0, newMaxVal], nRow - 1, 1);
        while newSumColumn(iColumn) - sum(tempColumn) > maxVal || newSumColumn(iColumn) - sum(tempColumn) < 0
            tempColumn = randi([0, newMaxVal], nRow - 1, 1);
        end
        outputMatrix(:, iColumn) = outputMatrix(:, iColumn) + [tempColumn; newSumColumn(iColumn) - sum(tempColumn)];
    end % for iColumn
else
    for iColumn = 1:nColumn
        nElement = newSumColumn(iColumn) + nRow;
        idPartition = [0, sort(randperm(nElement - 1, nRow - 1)), nElement];
        while any(diff(idPartition) - 1 > (newMaxVal))
            idPartition = [0, sort(randperm(nElement - 1, nRow - 1)), nElement];
        end
        outputMatrix(:, iColumn) = outputMatrix(:, iColumn) + reshape(diff(idPartition) - 1, [], 1);
    end
end
end % generaterandompartitionmatrixcolumn