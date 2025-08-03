"""
    QSMPhantom

A Julia module for Quantitative Susceptibility Mapping (QSM) phantom
generation and simulation including tissue properties, MRI signal models,
and data saving utilities.
"""
module QSMPhantom

using CSV, DataFrames, JLD2, ImageFiltering, FFTW, MeshGrid, Statistics, NIfTI

# Include tissue properties functions
include("tissue_properties/read_csv_model.jl")
include("tissue_properties/local_susceptibility_per_tissue.jl")
include("tissue_properties/local_susceptibility_exclude.jl")
include("tissue_properties/background_susceptibility.jl")
include("tissue_properties/r2star_per_tissue.jl")
include("tissue_properties/include_lobe.jl")
include("tissue_properties/partial_volume_effect.jl")

# Include MRI signal functions
include("mri_signal/anatomical_padding.jl")
include("mri_signal/anatomical_unpadding.jl")
include("mri_signal/dipole_kernel.jl")
include("mri_signal/dipole_convolution.jl")
include("mri_signal/demean.jl")
include("mri_signal/fat_peaks.jl")
include("mri_signal/simulate_acquisition.jl")
include("mri_signal/itok.jl")
include("mri_signal/ktoi.jl")

# Include data saving utilities
include("save_data/save_nii.jl")

# Export all functions in one statement for clarity
export read_csv_model, local_susceptibility_per_tissue, local_susceptibility_exclude,
       background_susceptibility, r2star_per_tissue, include_lobe, partial_volume_effect,
       anatomical_padding, anatomical_unpadding, dipole_kernel, dipole_convolution,
       demean, simulate_acquisition, fat_peaks, itok, ktoi, save_nii

end
