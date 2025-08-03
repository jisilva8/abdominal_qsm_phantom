"""
    dipole_convolution(chi_synth::Array{Float64,3}) -> Array{Float64,3}

Performs the convolution between the dipole kernel and the synthetic susceptibility map,
using the discrete kernel model from Milovic et al. (2017).

# Arguments
- `chi_synth::Array{Float64,3}`: 3D synthetic susceptibility map.

# Returns
- `phi_synth::Array{Float64,3}`: 3D fieldmap resulting from the convolution with the dipole kernel.

# Example
    phi = dipole_convolution(chi_synth)

# Notes
- Uses fixed voxel size (2, 2, 2) mm for kernel generation.
- Kernel model is set to discrete (Milovic et al. 2017).

Last edited by Silva J on 2025.08.03.
"""
function dipole_convolution(chi_synth::Array{Float64, 3})
    #Generate the dipole kernel
    matrix_size = size(chi_synth)
    kernel      = dipole_kernel(matrix_size, (2,2,2), 1)
    
    #Perform the convolution
    phi_synth   = real(ifft(kernel .* fft(chi_synth)))

    return phi_synth
end

