"""
    ktoi(data::Array; dims=nothing) -> Array

Computes the inverse Fourier transform of the input array `data` with optional
dimension selection. Applies `fftshift` and `ifftshift` before and after
transform to center the zero-frequency component.

# Arguments
- `data::Array`: Input array in k-space (frequency domain).
- `dims::Union{Nothing, AbstractVector{Int}}`: Optional vector of dimensions to apply the inverse FFT. Defaults to `nothing` (all dimensions).

# Returns
- `Array`: The inverse Fourier transformed array, shifted to keep zero-frequency 
component centered.

# Example
    img = ktoi(kspace_data)
    img2 = ktoi(kspace_data, dims=[1,2])

# Notes
- If `dims` is provided, the inverse FFT is applied sequentially along those 
dimensions.
- Uses `fftshift` and `ifftshift` to maintain spatial alignment.

Last edited by Silva J. (July 2025).
"""
function ktoi(data::Array, dims=nothing)
    if dims === nothing
        return fftshift(ifft(ifftshift(data)))
    else
        for d in dims
            data = fftshift(ifft(ifftshift(data, d), d), d)
        end
        return data
    end
end