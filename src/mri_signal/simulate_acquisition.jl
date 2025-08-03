"""
    simulate_acquisition(phi::Array{Float64,3}, R2s::Array{Float64,3}, TE_ms::Vector{Float64}, mask::Array{Float64,3}, Npeaks::Int=6, rBW::Float64=2000.0) -> Tuple{Array{ComplexF64,4}, Array{ComplexF64,4}, Array{ComplexF64,4}, Array{Float64,3}}

Simulates a multi-echo gradient echo (mGRE) acquisition using user-defined parameters. The simulation includes multipeak fat models, type-1 chemical shift artifact displacement, and variable SNR proportional to the receiver bandwidth.

# Arguments
- `phi::Array{Float64,3}`: 3D synthetic local or total fieldmap (ppm).
- `R2s::Array{Float64,3}`: 3D synthetic R2* map (Hz).
- `TE_ms::Vector{Float64}`: Echo times to be simulated (in milliseconds).
- `mask::Array{Float64,3}`: 3D binary mask of the region of interest (ROI).
- `Npeaks::Int=6`: Number of peaks considered in the fat model (default: 6). See 'fat_peaks.jl' for available peaks.
- `rBW::Float64=2000`: Receiver bandwidth in Hz/px (default: 2000).

# Returns
- `signal::Array{ComplexF64,4}`: Complex 4D matrix (Nx, Ny, Nz, Nte) with combined water and fat signal.
- `water_signal::Array{ComplexF64,4}`: Complex 4D matrix (Nx, Ny, Nz, Nte) with water-only signal.
- `fat_signal::Array{ComplexF64,4}`: Complex 4D matrix (Nx, Ny, Nz, Nte) with fat-only signal.
- `pdff::Array{Float64,3}`: Ground truth Proton Density Fat Fraction (PDFF).

# Example
    S, Sw, Sf, pdff = simulate_acquisition(phi, R2s, [2.3, 4.6, 6.9], mask, 3, 400.0)

# Notes
- The function uses preloaded experimental data for water, fat, and initial phase.
- Chemical shift displacement artifacts are simulated using a type-1 model.
- SNR scales inversely with the square root of receiver bandwidth.
- Noise is added based on the SNR model assuming optimal SNR=1000 at rBW=50 Hz/px.

Last edited by Silva J on 2025.08.03.
"""
function simulate_acquisition(phi::Array{Float64,3}, R2s::Array{Float64,3}, TE_ms::Vector{Float64},
     mask::Array{Float64,3}, Npeaks::Int=6, rBW::Float64=2000.0)

    #Load water, fat, and initial phase from experimental data
    water_file  = joinpath(@__DIR__, "..", "models", "common", "wtrexp.jld2")
    fat_file    = joinpath(@__DIR__, "..", "models", "common", "fatexp.jld2")
    phs0_file   = joinpath(@__DIR__, "..", "models", "common", "phszero.jld2")
    water_mag   =  jldopen(water_file, "r") do f; f["wtrexp"]; end
    fat_mag     =  jldopen(fat_file, "r") do f; f["fatexp"]; end
    phs_zero    =  jldopen(phs0_file, "r") do f; f["phszero"]; end
    
    #Parameters, constants, masking and simplified expressions
    B0          = 3;                 
    gyro        = 2 * pi * 42.58
    Nte         = length(TE_ms)
    TE_s        = TE_ms / 1000
    water_mag   = water_mag .* mask
    fat_mag     = fat_mag .* mask
    phs_zero    = phs_zero .* mask
    phi         = phi .* mask
    R2s         = R2s .* mask
    Nx,Ny,Nz    = size(phi)

    #Different initial phase for water and fat, based on pdff and pdwf
    pdwf = abs.(water_mag) ./ (abs.(water_mag) .+ abs.(fat_mag))
    pdff = abs.(fat_mag) ./ (abs.(water_mag) .+ abs.(fat_mag))
    pdwf[isnan.(pdwf)] .= 0
    pdff[isnan.(pdff)] .= 0
    phs_zero_water      = phs_zero .* pdwf
    phs_zero_fat        = phs_zero .* pdff

    #Compute fat signal from multipeak model
    da, df      = fat_peaks(Npeaks)
    signal_df   = da' .* exp.(1im * B0 * gyro * (TE_s * df'))
    
    #Define signal elements as 5D arrays with format (Nx,Ny,Nz,Nte,Npeaks)
    water_mag       = reshape(water_mag, Nx, Ny, Nz, 1, 1)
    fat_mag         = reshape(fat_mag, Nx, Ny, Nz, 1, 1)
    phs_zero_water  = reshape(phs_zero_water, Nx, Ny, Nz, 1, 1)
    phs_zero_fat    = reshape(phs_zero_water, Nx, Ny, Nz, 1, 1)
    pdwf            = reshape(pdwf, Nx, Ny, Nz, 1, 1)
    pdff            = reshape(pdff, Nx, Ny, Nz, 1, 1)
    signal_df       = reshape(signal_df, 1, 1, 1, Nte, Npeaks)
    R2s_decay       = R2s .* reshape(TE_s, 1, 1, 1, :, 1)
    phs_B0          = (B0 * gyro) .* phi .* reshape(TE_s, 1, 1, 1, :,1)

    #Compute fat displacement factor (Type 1 Chemical shift artifact)
    xx      = range(-1, stop=1 - 2/Nx, length=Nx)
    yy      = range(-1, stop=1 - 2/Ny, length=Ny)
    zz      = range(-1, stop=1 - 2/Nz, length=Nz)
    kspace, _, _ = meshgrid(yy, xx, zz)
    shift_factor = exp.(-1im * (B0 * gyro / rBW) .* kspace .* reshape(df, 1,1,1,1,:))

    #Compute water and fat signal
    water_signal    = pdwf .* water_mag .* exp.(1im .* (phs_B0 .+ phs_zero_water) - R2s_decay)
    fat_signal      = pdff .* fat_mag .* signal_df .* exp.(1im .* (phs_B0 .+ phs_zero_fat) - R2s_decay)
    
    #Apply Chemical shift displacement and combine the peaks
    fat_signal      = ktoi(itok(fat_signal, [1,2,3]) .* shift_factor, [1,2,3])
    fat_signal      = sum(fat_signal, dims=5)

    #Drop empty dimensions and combine signals}
    pdwf            = dropdims(pdwf, dims=(4,5))
    pdff            = dropdims(pdff, dims=(4,5))
    water_signal    = dropdims(water_signal, dims=5)
    fat_signal      = dropdims(fat_signal, dims=5) 
    
    #Add noise using a rBW-dependant SNR model
    #An optimal SNR=1000 at rBW=50 Hz/px is assumed
    base_SNR = 1000;
    base_rBW = 50;
    SNR   = base_SNR * sqrt(base_rBW / rBW)
    noise = (1 / SNR) * (randn(ComplexF64, size(phi)...))

    #Normalize
    signal          = water_signal + fat_signal
    water_signal    = water_signal ./ maximum(abs.(signal)) .+ noise
    fat_signal      = fat_signal ./ maximum(abs.(signal)) .+ noise
    signal          = water_signal + fat_signal

    return signal, water_signal, fat_signal, pdff
end