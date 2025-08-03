"""
    local_susceptibility_per_tissue(chi_model::Matrix{Float64}) 
        -> Tuple{Array{Float64,4}, Array{Float64,3}}

Assigns magnetic susceptibility values to each tissue of the phantom using 
realistic textures derived from anatomical priors. The resulting data is stored 
in a 4D array, with each tissue occupying a separate slice in the 4th dimension.

# Arguments
- `chi_model::Matrix{Float64}`: A (N × M) matrix of tissue model parameters, where 
  `N` is the number of tissues and `M` the number of model coefficients.

# Returns
- `chi_synth::Array{Float64,4}`: A 4D array of size (X, Y, Z, N) containing synthetic 
  susceptibility values for each of the `N` tissues.
- `chi_mask::Array{Float64,3}`: A 3D binary mask indicating phantom tissues (`1`) 
  versus background air (`0`).

# Example
    chi_synth, chi_mask = local_susceptibility_per_tissue(chi_model)

# Notes
- Textures are modulated using precomputed anatomical maps: `r2snorm`, `wtrnorm`, 
  and `fatnorm`.
- The final susceptibility values are generated as a linear combination of 
  texture maps and tissue-specific coefficients.
- The first column of `chi_model` provides the baseline susceptibility for each tissue.

Last edited by Silva J on 2025.08.03.
"""
function local_susceptibility_per_tissue(chi_model::Matrix{Float64})
    
    #Load tissue properties (labels, and normalized maps)
    
    tlabels_file = joinpath(@__DIR__, "..", "models", "common", "tissue_masks.jld2")
    r2snorm_file = joinpath(@__DIR__, "..", "models", "common", "r2snorm.jld2")
    wtrnorm_file = joinpath(@__DIR__, "..", "models", "common", "wtrnorm.jld2")
    fatnorm_file = joinpath(@__DIR__, "..", "models", "common", "fatnorm.jld2")
    
    tissue_masks = jldopen(tlabels_file, "r") do f; f["tissue_masks"]; end
    r2snorm      = jldopen(r2snorm_file, "r") do f; f["r2snorm"]; end
    wtrnorm      = jldopen(wtrnorm_file, "r") do f; f["wtrnorm"]; end
    fatnorm      = jldopen(fatnorm_file, "r") do f; f["fatnorm"]; end
    
    #Relevant sizes
    phantom_size = size(tissue_masks)
    n_tissue     = maximum(tissue_masks)

    #Assign susceptibility values to each tissue
    chi_synth = zeros((phantom_size...,n_tissue))
    for tissue in 1:n_tissue
        #Mean susceptibility value 
        chi_mean = chi_model[tissue,1]
        
        #Texture modulation
        r2s_mod = chi_model[tissue, 5] .* (r2snorm .- chi_model[tissue, 2]);
        wtr_mod = chi_model[tissue, 6] .* (wtrnorm .- chi_model[tissue, 3]);
        fat_mod = chi_model[tissue, 7] .* (fatnorm .- chi_model[tissue, 4]);

        #Susceptibility map for tissue
        chi_synth[:,:,:,tissue] .= chi_mean .+ r2s_mod .+ wtr_mod .+ fat_mod;
    end

    #Phantom mask with external air excluded
    chi_mask = Float64.(tissue_masks .> 0)
    
    return chi_synth, chi_mask
end