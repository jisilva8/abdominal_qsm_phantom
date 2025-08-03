"""
    r2star_per_tissue(r2star_model::Matrix{Float64}) -> Tuple{Array{Float64,4}, Array{Float64,3}}

Generates R2* maps with realistic textures for each tissue in the phantom.  
Each modified tissue is stored separately along the 4th dimension of a 4D array.

# Arguments
- `r2star_model::Matrix{Float64}`: A numeric matrix of size (N × M), where `N` is the number of tissues, and `M` contains parameters 
  for whether and how each tissue is modified.

# Returns
- `r2star_synth::Array{Float64,4}`: 4D array of shape (X, Y, Z, N), where each slice along the 4th dimension contains R2* values 
  for one modified tissue. The first slice holds the unmodified base.
- `r2star_mask::Array{Float64,3}`: 3D array labeling each voxel by its modified tissue ID. Used for partial volume modeling.

# Example
    r2star_synth, r2star_mask = r2star_per_tissue(r2star_model)

# Notes
- Only tissues marked with `1` in the first column of `r2star_model` are modified.
- R2* textures are generated using normalized templates (r2snorm) and modulation parameters.
- Tissue masks and base R2* maps are loaded from JLD2 files in the `models/common/` directory.

Last edited by Silva J on 2025.08.03.
"""
function r2star_per_tissue(r2star_model::Matrix{Float64})
    
    #Load tissue properties (labels, and normalized maps)
    tlabels_file = joinpath(@__DIR__, "..", "models", "common", "tissue_masks.jld2")
    r2snorm_file = joinpath(@__DIR__, "..", "models", "common", "r2snorm.jld2")
    wtrnorm_file = joinpath(@__DIR__, "..", "models", "common", "wtrnorm.jld2")
    fatnorm_file = joinpath(@__DIR__, "..", "models", "common", "fatnorm.jld2")
    r2sreal_file = joinpath(@__DIR__, "..", "models", "common", "r2sreal.jld2")

    tissue_masks = jldopen(tlabels_file, "r") do f; f["tissue_masks"]; end
    r2sreal      = jldopen(r2sreal_file, "r") do f; f["r2sreal"]; end
    r2snorm      = jldopen(r2snorm_file, "r") do f; f["r2snorm"]; end
    wtrnorm      = jldopen(wtrnorm_file, "r") do f; f["wtrnorm"]; end
    fatnorm      = jldopen(fatnorm_file, "r") do f; f["fatnorm"]; end
    
    #Relevant sizes
    phantom_size = size(tissue_masks)
    n_tissue     = maximum(tissue_masks)

    #Generate initial R2* map and R2* mask
    r2star_synth = reshape(r2sreal, size(r2sreal)..., 1)
    r2star_mask  = Float64.(tissue_masks .> 0)

    #Keep track of modified tissues for future partial volumes
    # 0 = air region
    # 1 = unmodified tissue
    # N = (N-1)th  modified tissue
    mod_id  = 1

    #Assign R2* values and R2* mask
    for tissue in 1:n_tissue
        #Modify only selected tissues
        if r2star_model[tissue,1] == 1
            #Create mask for partial volume effects
            mod_id      += 1
            tmp_mask     = Float64.(tissue_masks .== tissue) 
            r2star_mask  = r2star_mask .* (1 .- tmp_mask) + tmp_mask .* mod_id
            
            #Modify selected value from current R2* map
            r2star_mean  = r2star_model[tissue, 2] 
            r2star_mod   = r2star_model[tissue, 4] .* (r2snorm .- r2star_model[tissue, 3]);
            r2star_synth = cat(r2star_synth, r2star_mean .+ r2star_mod; dims=4)
        end
    end
    
    return r2star_synth, r2star_mask
end