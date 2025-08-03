"""
    background_susceptibility(chi_model::Matrix{Float64}, chi_air::Float64, exclusions::Vector{Int64}=Int64[]) -> Array{Float64,3}

Generates a 3D synthetic susceptibility map of background sources.

# Arguments
- `chi_model::Matrix{Float64}`: Matrix of size (N x M), where N is the number of tissues and M is the number of parameters of the susceptibility model.
- `chi_air::Float64`: Susceptibility value of external air, in ppm.
- `exclusions::Vector{Int64}`: List of tissue IDs to be considered as background sources (e.g. 12 = large intestine).

# Returns
- `chi_synth::Array{Float64,3}`: 3D susceptibility map including background sources.

# Example
    chi_synth_background = background_susceptibility(chi_model, chi_air, exclusions)

# Notes
- Background sources include internal air cavities like lung, stomach, and colon.
- External air susceptibility is included as a separate label with tissue ID 0.
- The function sums contributions from all excluded tissues and external air.

Last edited by Silva J on 2025.08.03.
"""
function background_susceptibility(chi_model::Matrix{Float64}, chi_air::Float64, exclusions::Vector{Int64} = Int64[])

    #Load tissue labels
    tlabels_file = joinpath(@__DIR__, "..", "models", "common", "tissue_masks.jld2")
    tissue_masks = jldopen(tlabels_file, "r") do f; f["tissue_masks"]; end

    #Relevant sizes
    phantom_size = size(tissue_masks)
    n_tissue     = maximum(tissue_masks)

    #Generate susceptibility map of internal air cavities (e.g., lung, stomach, colon)
    chi_synth = zeros((phantom_size...,n_tissue))
    for tissue in exclusions
        chi_synth[:,:,:,tissue] .= Float64.(tissue_masks .== tissue) .* chi_model[tissue,1]
    end

    #Include external air
    chi_synth = cat(chi_synth, Float64.(tissue_masks .== 0) .* chi_air; dims=4)    
    chi_synth = sum(chi_synth, dims=4)
    chi_synth = dropdims(chi_synth, dims = 4)

    return chi_synth
end