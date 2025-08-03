"""
    partial_volume_effect(tp_synth::Array{Float64,4}, sigma::Float64, exclusions::Vector{Int}=Int[]) -> Array{Float64,3}

Applies partial volume effects to a synthetic tissue property map by smoothing 
the interfaces between tissue types using Gaussian filtering.

# Arguments
- `tp_synth::Array{Float64,4}`: 4D array of shape (X, Y, Z, N) containing synthetic values 
  for N different tissues. Each tissue occupies one slice in the 4th dimension.
- `sigma::Float64`: Standard deviation of the Gaussian kernel used for filtering.
- `exclusions::Vector{Int}` (optional): List of tissue IDs to exclude from the filtering process.

# Returns
- `tp_synth_pv::Array{Float64,3}`: 3D array representing the synthetic property map with 
  partial volume effects applied to non-excluded tissues.

# Example
    chi_synth_pv = partial_volume_effect(chi_synth, sigma, [12, 23])

# Notes
- Only tissues not listed in `exclusions` are smoothed and blended.
- The function assumes that each tissue occupies its own 4D slice in `tp_synth`.
- Gaussian filtering is applied over tissue borders to simulate realistic transitions.

Last edited by Silva J on 2025.08.03.
"""
function partial_volume_effect(tp_synth::Array{Float64, 4}, sigma::Float64; 
    exclusions::Vector{Int64} = Int64[], tissue_masks::Union{Nothing, Array{<:Number, 3}} = nothing)
    

    #If no labels are provided, load default tissue labels
    if isnothing(tissue_masks)
        tlabels_file = joinpath(@__DIR__, "..", "models", "common", "tissue_masks.jld2")   
        tissue_masks = jldopen(tlabels_file, "r") do f; f["tissue_masks"]; end
    end

    #Relevant sizes
    phantom_size = size(tissue_masks)
    n_tissue     = size(tp_synth, 4)

    #List of excluded tissues
    exclusion_flags =  zeros(Int, n_tissue)
    exclusion_flags[exclusions] .= 1;

    #Generate probability masks for each tissue
    prob_mask   = zeros((phantom_size...,n_tissue))
    kernel      = KernelFactors.gaussian((sigma,sigma,sigma))

    for tissue in 1:n_tissue
        mask = Float64.(tissue_masks .== tissue)

        #For excluded tissues: Normal mask
        if exclusion_flags[tissue] == 1
            prob_mask[:,:,:,tissue] = mask
        
        #For non-excluded tissues: Gausian edge smoothed mask
        else
            prob_mask[:,:,:,tissue] = imfilter(mask, kernel)
        end
    end

    #Normalize probabilities considering the prob_mask overlaps
    prob_mask .= prob_mask ./ sum(prob_mask, dims=4)
    prob_mask .= ifelse.(isnan.(prob_mask), 0.0, prob_mask)

    #Combine all the tissues into a 3D tp_synth array with partial volumes
    tp_synth_pv = sum(prob_mask .* tp_synth, dims = 4)
    tp_synth_pv = dropdims(tp_synth_pv, dims = 4)

    return tp_synth_pv
end