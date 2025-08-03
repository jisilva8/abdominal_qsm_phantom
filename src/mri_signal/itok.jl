"""
    itok(data::Array; dims=nothing) -> Array

Computes the forward Fourier transform of the input array `data` with optional
dimension selection. Applies `ifftshift` and `fftshift` before and after
transform to center the zero-frequency component.

# Arguments
- `data::Array`: Input array in image (spatial) domain.
- `dims::Union{Nothing, AbstractVector{Int}}`: Optional vector of dimensions to apply the FFT. Defaults to `nothing` (all dimensions).

# Returns
- `Array`: The Fourier transformed array, shifted to keep zero-frequency component centered.

# Example
    kspace = itok(image_data)
    kspace2 = itok(image_data, dims=[1,2])

# Notes
- If `dims` is provided, the FFT is applied sequentially along those dimensions.
- Uses `fftshift` and `ifftshift` to maintain frequency alignment.

Last edited by Silva J. (July 2025).
"""
function itok(data::Array, dims=nothing)
    if dims === nothing
        return fftshift(fft(ifftshift(data)))
    else
        for d in dims
            data = fftshift(fft(ifftshift(data, d), d), d)
        end
        return data
    end
end