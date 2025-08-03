"""
    anatomical_padding(chi_synth::Array{Float64,3}) -> Array{Float64,3}

Triples the size of the `chi_synth` array. Padding along the x and y dimensions
is done by filling with the susceptibility value of the surrounding air.
Padding along the z dimension is performed by appending inverted copies of the
array to ensure anatomically consistent continuity.

# Arguments
- `chi_synth::Array{Float64,3}`: (x,y,z) array with synthetic susceptibility values.

# Returns
- `chi_pad::Array{Float64,3}`: 3D array padded along the three dimensions.

# Example
    chi_pad = anatomical_padding(chi_synth)

Last edited by Silva J. (July 2025).
"""
function anatomical_padding(chi_synth::Array{Float64, 3}) 
    
    #Relevant sizes
    phantom_size = size(chi_synth)
    
    #Padding along X and Y, instead of zeros, we will padding
    #with the air susceptibility value
    chi_pad_xy = padarray(chi_synth, Fill(chi_synth[1,1,1], (192,192,0)))

    #Generate reversed copies and append them along Z axis
    chi_flip = chi_pad_xy[:, :, end:-1:1]
    chi_pad = cat(chi_flip, chi_pad_xy, chi_flip; dims=3)
    
    return chi_pad
end