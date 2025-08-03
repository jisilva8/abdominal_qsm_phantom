"""
    fat_peaks(Npeaks::Int) -> Tuple{Vector{Float64}, Vector{Float64}}

Provides the fat model parameters for a given number of fat peaks.

# Arguments
- `Npeaks::Int`: Number of fat peaks in the model.

# Returns
- `fat_ramps::Vector{Float64}`: Relative amplitudes of the fat peaks.
- `fat_freqs::Vector{Float64}`: Relative frequencies of the fat peaks (in ppm), 
  with respect to water.

# Example
    ramps, freqs = fat_peaks(6)

# Notes
- Supported Npeaks values: 1, 3, 4, 5, 6, 7, 9.
- If an unsupported number is given, the function throws an error.

Last edited by Silva J on 2025.08.03.
"""
function fat_peaks(Npeaks::Int)
    fat_ramps = []
    fat_freqs = []

    if Npeaks == 1
        fat_ramps = [1.0]
        fat_freqs = [-3.4]
    elseif Npeaks == 3
        fat_ramps = [0.080, 0.170, 0.750]
        fat_freqs = [0.73, -2.49, -3.29]
    elseif Npeaks == 4
        fat_ramps = [0.080, 0.150, 0.720, 0.040]
        fat_freqs = [0.73, -2.49, -3.27, -3.68]
    elseif Npeaks == 5
        fat_ramps = [0.080, 0.050, 0.100, 0.720, 0.040]
        fat_freqs = [0.73, -2.35, -2.54, -3.27, -3.68]
    elseif Npeaks == 6
        fat_ramps = [0.088, 0.700, 0.120, 0.006, 0.039, 0.047]
        fat_freqs = [-3.80, -3.40, -2.60, -1.95, -0.50, 0.60]
    elseif Npeaks == 7
        fat_ramps = [0.042, 0.015, 0.066, 0.096, 0.071, 0.627, 0.083]
        fat_freqs = [0.61, -1.93, -2.45, -2.67, -3.11, -3.4, -3.8]
    elseif Npeaks == 9
        fat_ramps = [0.037, 0.010, 0.039, 0.006, 0.058, 0.062, 0.058, 0.642, 0.088]
        fat_freqs = [0.59, 0.49, -0.5, -1.95, -2.46, -2.68, -3.1, -3.4, -3.8]
    else
        error("No valid value for Npeaks")
    end

    return fat_ramps, fat_freqs
end