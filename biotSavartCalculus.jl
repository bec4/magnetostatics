using ForwardDiff

function pathVec(pointArray::Array{Float64,2})
    #assume pointArray is an 2d array 3 X N
    #assume the path vector associated with the end of the path is a zero vector.
    return hcat(pointArray[:,2:end]-pointArray[:,1:end-1],[0.0, 0.0, 0.0])
end


function normCubed(x)
    #calculates |x|^3
    return sqrt.(sum(x.^2, dims = 1)).^3
end

function crossArray(x, y)
    #assume input arrays are 2d array 3 X N
    #returns 3XN array
    #return hcat(x[2,:].*y[3,:]-y[2,:].*x[3,:], x[3,:].*y[1,:]-y[3,:].*x[1,:],x[1,:].*y[2,:] - y[1,:].*x[2,:])'
    xX = @view x[1,:]
    xY = @view x[2,:]
    xZ = @view x[3,:]
    yX = @view y[1,:]
    yY = @view y[2,:]
    yZ = @view y[3,:]
    return hcat(xY.*yZ-yY.*xZ, xZ.*yX-yZ.*xX,xX.*yY - yX.*xY)'
end


function biasField(probePoint, pathElement, dlvec)
    #computes dl × r./normCubed(rvec)
    rvec = pathElement.-probePoint
    return sum(crossArray(dlvec,rvec)./normCubed(rvec), dims = 2)
end

function fieldNorm(x, pathElement, dlvec)
   return sqrt(sum(biasField(x,pathElement,dlvec).^2))
end

function fieldZComponent(x, pathElement, dlvec)
    return biasField(x, pathElement, dlvec)[3]
end

function fieldNormplusBias(x, pathElement, dlvec, bias = [0, 0, 0])
   return sqrt(sum( ( biasField(x,pathElement,dlvec) + bias  ).^2))
end