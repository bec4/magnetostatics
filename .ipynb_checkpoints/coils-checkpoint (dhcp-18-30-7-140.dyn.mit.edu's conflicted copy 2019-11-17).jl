# Function that returns the path coordinates of a regularly wound coil. First we spiral up, the number of windings equal to
# `pitches,` then we spiral down, repeating this process `nlayer` times. 
function normalCoil(startingAngle, pitches, nlayers, bucketWallRadius, sign, wt=3.5, angRes=0.003, zpos = 0)
    # Preallocate
    coil = Array{Float64,2}(undef, 0, 3);
    
    endingAngle = 2*pi*pitches+startingAngle;
    ang = collect(range(startingAngle, endingAngle, length = Int(floor(2*pi*pitches/angRes))));
    
    # Angles associated with top/bottom layers of each spiral
    extremaAng = range(startingAngle + angRes, startingAngle + 2*pi, step = angRes);

    for layer=1:nlayers
        radius = bucketWallRadius + (layer-1)*wt + wt/2;
       
        # Winding up
        if mod(layer,2)==1
            tempcoil=radius.*hcat(cos.(ang), sign.*sin.(ang), (ang .- startingAngle)./(2*pi*radius).*wt.+zpos./radius);
            
            # Don't add extra top winding if this is the last layer
            if layer ≠ nlayers
                topWinding = hcat( (radius .+ (extremaAng .- startingAngle) * wt / (2*pi)) .* cos.(extremaAng), 
                    (radius .+ (extremaAng .- startingAngle) * wt /(2*pi)) .* sin.(extremaAng).*sign,
                    pitches*wt .* ones(length(extremaAng)) .+ zpos);
                tempcoil = vcat(tempcoil, topWinding);
            end
        
        # Winding down
        else
            tempcoil=radius.*hcat(cos.(ang), sign.*sin.(ang), (endingAngle.-ang)./(2*pi*radius).*wt.+zpos./radius);
            
            # Don't add extra bottom winding if this is the last layer
            if layer ≠ nlayers
                bottomWinding = hcat( (radius .+ extremaAng * wt / (2*pi)) .* cos.(extremaAng),
                    (radius .+ extremaAng * wt /(2*pi)) .* sin.(extremaAng).*sign,
                    zpos .* ones(length(extremaAng)));
                tempcoil = vcat(tempcoil, bottomWinding);
            end
        end
 
        coil = vcat(coil, tempcoil);
        
    end
    return coil
end

# Function that returns the path coordinates for a bilayer, counterwound coil, which minimizes the curvature at the
# center in a Helmholtz configuration. By definition, there are two layers, with opposite helicity.
function coilStack(startingAngle, turns, bucketWallRadius, sign, wt=3.5, angRes=0.003, zpos = 0)
    ang = collect(range(startingAngle, 2*pi*turns+startingAngle, length = Int(floor(2*pi*turns/angRes))))
    ang2 = collect(range(2*pi*turns+startingAngle,2*pi*2*turns+startingAngle, length = Int(floor(2*pi*turns/angRes))))
    
    #generate the top spiral of the double layer, inward
    radius = bucketWallRadius + turns*wt+wt/2; # outer radius, wt/2 is the center of the coil
    coreTop = (radius.-ang/(2*pi) * wt).*hcat(cos.(ang), sign.*sin.(ang), zeros(size(ang)));
    coreTop[:,3] = zpos .+ fill(wt,length(ang),1)

    #generate the bottom spiral of the double layer, outward
    innerRadius = radius-ang[end]/(2*pi)*wt;
    transitionLength = Int(floor(0.2/angRes)); #transition part is about 11 degrees
    coreBottom = (innerRadius.+(ang2.-ang[end])/(2*pi)*wt).*hcat(cos.(ang),sign.*sin.(ang), zeros(size(ang)))

    #generate the transition part, at the transition part, the coil form a slope in Z component
    coreBottom[1:transitionLength,3] = ( zpos + wt ) .- wt/transitionLength *(1:1.0:transitionLength)

    coreBottom[transitionLength + 1 : end, 3] = zpos * ones( length(coreBottom[transitionLength + 1 : end, 3]) )
        
    return vcat(coreTop, coreBottom)
end