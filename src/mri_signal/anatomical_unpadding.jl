"""
    anatomical_unpadding(chi_pad::Array{Float64,3}) -> Array{Float64,3}

Undoes the padding process from `anatomical_padding(chi_synth)`, recovering the original dimensions of `chi_synth`.

# Arguments
- `chi_pad::Array{Float64,3}`: 3D array padded along the three dimensions.

# Returns
- `chi_synth::Array{Float64,3}`: Unpadded array with its original dimensions.

# Example
    chi_synth = anatomical_unpadding(chi_pad)

Last edited by Silva J. (July 2025).
"""
function anatomical_unpadding(chi_pad::Array{Float64,3})
    
    #Rescue the original size before padding
    chi_unpad = chi_pad[(end÷3 + 1):(2end÷3), (end÷3 + 1):(2end÷3), (end÷3 + 1):(2end÷3)]

    return chi_unpad
end