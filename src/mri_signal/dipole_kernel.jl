"""
    dipole_kernel(matrix_size::Tuple{Int,Int,Int}, voxel_size::Tuple{Int,Int,Int}, kernel_model::Int) -> Array{Float64,3}

Calculates the dipole kernel used in QSM (single orientation) that models the susceptibility-to-field convolution.
Use this function if the main magnetic field direction is along the z axis.
Includes both continuous and discrete kernel approximations.

# Arguments
- `matrix_size::Tuple{Int,Int,Int}`: Dimensions of the 3D array.
- `voxel_size::Tuple{Int,Int,Int}`: Spatial resolution of each voxel.
- `kernel_model::Int`: Dipole model selector:
    - 0: Continuous kernel (Salomir et al. 2003).
    - 1: Discrete kernel (Milovic et al. 2017).

# Returns
- `kernel::Array{Float64,3}`: 3D array with the dipole kernel in frequency space.

# Example
    kernel = dipole_kernel((64,64,64), (1,1,1), 1)

# Notes
- The function shifts the kernel center to [1,1,1] in frequency domain.
- Kernel center value is set to zero to avoid division by zero artifacts.

Created by Milovic C on 30.03.2017.  
Last edited by Silva J on 2025.08.03.
"""
function dipole_kernel(matrix_size::Tuple{Int64, Int64, Int64}, voxel_size::Tuple{Int64, Int64, Int64},
    kernel_model::Int64)

    Nx,Ny,Nz = matrix_size

    #Continous kernel model (Salomir et al. 2003)
    if kernel_model == 0
        #Generate grid for kspace
        ky_range = -floor(Int, Nx/2):ceil(Int, Nx/2)-1
        kx_range = -floor(Int, Ny/2):ceil(Int, Ny/2)-1
        kz_range = -floor(Int, Nz/2):ceil(Int, Nz/2)-1
        ky, kx, kz = meshgrid(ky_range, kx_range, kz_range)

        #K-space normalization and scaling
        kx = (Float32.(kx) ./ maximum(abs.(vec(Float32.(kx))))) / voxel_size[1]
        ky = (Float32.(ky) ./ maximum(abs.(vec(Float32.(ky))))) / voxel_size[2]
        kz = (Float32.(kz) ./ maximum(abs.(vec(Float32.(kz))))) / voxel_size[3]

        #Define the k² term
        k2 = kx.^2 + ky.^2 + kz.^2;
        k2[k2 .== 0] .= eps()

        #Kernel model
        kernel = 1/3 .- (kz .^ 2) ./ k2

    #Discrete kernel model (Milovic et al. 2017)
    else
        #Define kspace range and resolution
        cx, cy, cz = 1 .+ fld.(matrix_size, 2) 
        kx, ky, kz = (1:Nx) .- cx , (1:Ny) .- cy , (1:Nz) .- cz
        dx, dy, dz = 1 ./ (matrix_size .* voxel_size)


        #Set k-space
        kx, ky, kz = (kx, ky, kz) .* (dx,dy,dz)
        kx, ky, kz = (reshape(kx,:,1,1), reshape(ky,1,:,1), reshape(kz,1,1,:))
        kx, ky, kz = (repeat(kx,1,Ny,Nz), repeat(ky,Nx,1,Nz), repeat(kz,Nx,Ny,1))


        #Define the k² term
        k2 = -3 .+ cos.(2 * pi .* kx) .+ cos.(2 * pi .* ky) .+ cos.(2 * pi .* kz)
        k2[k2 .== 0] .= eps()

        #Kernel model
        kernel = 1/3 .- (-1 .+ cos.(2 * pi .* kz)) ./ k2

    end
    
    #Keep the center of the frequency domain at [1,1,1]
    kernel = ifftshift(kernel)
    kernel[1, 1, 1] = 0.0

    return kernel
end