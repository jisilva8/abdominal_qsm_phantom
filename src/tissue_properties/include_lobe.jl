"""
    include_lobe(tp_synth::Array{Float64,4}, lobe_mean::Float64=0.2, lobe_alpha::Float64=0.5; tissue_masks::Union{Nothing, Array{<:Number,3}}=nothing) -> Tuple{Array{Float64,4}, Array{Float64,3}}

Includes a synthetic lobe within the liver, with user-defined mean tissue property value and texture.

# Arguments
- `tp_synth::Array{Float64,4}`: 4D array of shape (X, Y, Z, N) containing synthetic values for N different tissues.
- `lobe_mean::Float64=0.2`: Mean value for the synthetic lobe.
- `lobe_alpha::Float64=0.5`: Texture modulation parameter.
- `tissue_masks::Union{Nothing, Array{<:Number,3}}` (optional): 3D labeled set of tissue masks. If not provided, default tissue labels are loaded.

# Returns
- `tp_synth_updated::Array{Float64,4}`: 4D array with shape (X, Y, Z, N+1) including the additional synthetic lobe.
- `tp_masks_updated::Array{Float64,3}`: 3D array with updated tissue masks including the new lobe label.

# Example
    tp_synth_updated, tp_masks_updated = include_lobe(tp_synth, 0.3, 0.7)

# Notes
- The function loads default tissue masks if `tissue_masks` is not provided.
- The new lobe is added as an additional tissue label with ID one greater than the maximum existing label.
- Texture modulation applies spatial variation based on a predefined texture mask.

Last edited by Silva J on 2025.08.03.
"""
function include_lobe(tp_synth::Array{Float64, 4}, lobe_mean::Float64 = 0.2, lobe_alpha::Float64 = 0.5;
    tissue_masks::Union{Nothing, Array{<:Number, 3}} = nothing)
    #Load lobe file
    lobe_file       = joinpath(@__DIR__, "..", "models", "common", "lobe.jld2") 
    lobe_mask       = jldopen(lobe_file, "r") do f; f["lobe_mask"]; end
    lobe_texture    = jldopen(lobe_file, "r") do f; f["lobe_texture"]; end

    #If no labels are provided, load default tissue labels
    if isnothing(tissue_masks)
        tlabels_file = joinpath(@__DIR__, "..", "models", "common", "tissue_masks.jld2")   
        tissue_masks = jldopen(tlabels_file, "r") do f; f["tissue_masks"]; end
    end

    #Include label for lobe
    tp_masks_updated = Float64.(tissue_masks) .* (1 .- lobe_mask)
    tp_masks_updated = tp_masks_updated + (lobe_mask .* (maximum(tissue_masks)+1))

    #Redefine susceptibility values for lobe within and outside the mask    
    lobe_tp = (lobe_mean .+ (lobe_alpha .* lobe_texture)) .* lobe_mask
    lobe_tp = lobe_tp + (1 .- lobe_mask) .* lobe_mean
    tp_synth_updated = cat(tp_synth, lobe_tp; dims=4)

    return tp_synth_updated, tp_masks_updated
end
