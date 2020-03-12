classdef LinearEncodingMatrix
%

properties (SetAccess = private, GetAccess = private)

    adjoint =   0;
end

properties (SetAccess = private, GetAccess = public)
    EMat   =   [];
    objectDims    =   [];
    Nenc      =   0;
    Ncomp      =  0;
end

methods
function res = LinearEncodingMatrix(EMat, dims)

    res.objectDims =   dims ;
    res.EMat =   EMat;
    res.Nenc =   size(EMat,1);
    res.Ncomp =   size(EMat,2);

end

function res = ctranspose(a)

    a.adjoint   = xor(a.adjoint,1);
    res         = a;

end

function res = mtimes(a,b)

    if a.adjoint
        % avoiding permute() by doing (x'A)' instead of A'x
        res = reshape(b, [], a.Nenc);
        res = res * a.EMat; %this will only work as long as VEMat is a real matrix
        res = reshape(res, a.objectDims(1), a.objectDims(2), a.objectDims(3), a.objectDims(4),[]);
    else
        % avoiding permute() by doing (x'A')' instead of Ax
        res = reshape(b, [], a.Ncomp);
        res = res * a.EMat'; %this will only work as long as VEMat is a real matrix
        res = reshape(res, a.objectDims(1), a.objectDims(2), a.objectDims(3), a.objectDims(4),[]);
    end

end


end
end
