"""
    demean(I::Array{Float64,3}, mask::Array{Float64,3}) -> Array{Float64,3}

Subtracts the mean value of the image `I` within a region of interest specified by `mask`.

# Arguments
- `I::Array{Float64,3}`: 3D image array to be demeaned.
- `mask::Array{Float64,3}`: 3D binary mask array defining the region of interest.

# Returns
- `I_demeaned::Array{Float64,3}`: Demeaned image with the mean value (within the mask) subtracted.

# Example
    I_demeaned = demean(I, mask)

# Notes
- The output image is zeroed outside the mask region.

Last edited by Silva J on 2025.08.03.
"""
function demean(I::Array{Float64, 3}, mask::Array{Float64, 3})

    #Remove mean value and return demeaned image
    I_demeaned = I .- mean(I[mask .== 1])
    I_demeaned = I_demeaned .* mask
    
    return I_demeaned

end