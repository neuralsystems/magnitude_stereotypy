function data = removeallnandims(data)
nDims = ndims(data);
function allnandims = findalldim(data, dims)
    if isempty(dims)
        allnandims = data;
    else
        allnandims = findalldim(all(data, dims(1)), dims(2:end));
    end
end
for iDim = 1:nDims
    subs_data = repelem({':'}, 1, nDims);
    subs_data{iDim} = ~findalldim(isnan(data), setdiff(1:nDims, iDim));
    data = subsref(data, substruct('()', subs_data));
end
end