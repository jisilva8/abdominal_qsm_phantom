"""
    local_susceptibility_exclude(
        chi_synth::Array{Float64,4}, 
        chi_mask::Array{Float64,3}, 
        exclusions::Vector{Int64}
    ) -> Tuple{Array{Float64,4}, Array{Float64,3}}

Removes selected tissues (considered as air) from a synthetic susceptibility 
phantom by zeroing out their values and updating the tissue mask accordingly.

# Arguments
- `chi_synth::Array{Float64,4}`: A 4D array of size (X, Y, Z, N) containing 
  synthetic susceptibility values for `N` tissues.
- `chi_mask::Array{Float64,3}`: A 3D binary mask indicating phantom tissues (`1`) 
  and background air (`0`).
- `exclusions::Vector{Int64}`: List of tissue IDs to exclude (e.g., `[12, 23]`).

# Returns
- `chi_synth::Array{Float64,4}`: Updated susceptibility map with excluded tissues set to `0`.
- `chi_mask::Array{Float64,3}`: Updated tissue mask with excluded tissues removed.

# Example
    chi_synth, chi_mask = local_susceptibility_exclude(chi_synth, chi_mask, [12, 23])

# Notes
- This function uses the original labeled tissue mask from `tissue_masks.jld2` to 
  identify which voxels to exclude.
- The susceptibility values of excluded tissues are set to `0` in all voxels.
- Excluded voxels are removed from the mask via element-wise multiplication.

Last edited by Silva J on 2025.08.03.
"""
function local_susceptibility_exclude(chi_synth::Array{Float64, 4}, chi_mask::Array{Float64, 3},
  exclusions::Vector{Int64})

    #Load tissue labels
    tlabels_file = joinpath(@__DIR__, "..", "models", "common", "tissue_masks.jld2") 
    tissue_masks = jldopen(tlabels_file, "r") do f; f["tissue_masks"]; end
    
    #Exclude tissues from list
    for tissue in exclusions        
        #Force excluded tissues' susceptibility to 0
        chi_synth[:,:,:,tissue] .= 0;

        #Update chi_mask
        chi_mask = chi_mask .* Float64.(1 .- (tissue_masks .== tissue))
    end

    return chi_synth, chi_mask
end