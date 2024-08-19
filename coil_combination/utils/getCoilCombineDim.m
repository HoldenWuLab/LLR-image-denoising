function order = getCoilCombineDim(dimStr)

if ~contains(dimStr,'c')
    error('It should contain coil dimension')
end

orderRef = 'xyznca';
Nstr = length(dimStr);
order = zeros(1,Nstr);

count = Nstr+1;
for i = 1:6
    strLoc = strfind(dimStr,orderRef(i));
    if ~isempty(strLoc)
        order(i) = strLoc;
    else
        order(i) = count;
        count = count + 1;
    end
end
