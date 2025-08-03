# =========================================================================================
# example03_synthetic_fieldmap.jl
# =========================================================================================
#
# Example script for generating synthetic local and total fieldmaps using QSMPhantom.
# These outputs can be used as ground truth for testing background field removal, 
# phase unwrapping, and dipole inversion algorithms.
#
# To run this script from the project root directory:
# - In Julia REPL:
#       include("examples/example03_synthetic_fieldmap.jl")
# - From the command line:
#       julia --project=. -e 'include("examples/example03_synthetic_fieldmap.jl")'
#
# Ensure the script is executed from the root folder to maintain correct relative paths.
#
# =========================================================================================
# Last edited by Silva J on 2025.08.03.
# =========================================================================================

# ====== Activate the project environment ======
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using NIfTI
using QSMPhantom

# ========= User-adjustable parameters ==========
output_folder       = joinpath(@__DIR__, "..", "results", "ex03_synthetic_fieldmap")
case_select         = 1             # 1 = healthy_subject; 2 = iron_overload; 3 = pathological_lobe
exclusions          = [12, 21, 23]  # Exclude heart (12), large intestine (21), and lungs (23)

# =================== MAIN CODE =================
# Select susceptibility model
cases                   = ["healthy_subject", "iron_overload", "pathological_lobe"]
chi_model_file          = joinpath(@__DIR__, "..", "src", "models", cases[case_select], "chi_model.csv")

# Generate local susceptibility
chi_model               = read_csv_model(chi_model_file)
chi_synth, chi_mask     = local_susceptibility_per_tissue(chi_model)
chi_synth, chi_mask     = local_susceptibility_exclude(chi_synth, chi_mask, exclusions)

# Generate background susceptibility
air_chi_ppm             = 9.2
chi_background          = background_susceptibility(chi_model, air_chi_ppm, exclusions)

# Add lobes (if selected) and apply partial volume effects
if case_select == 3
    chi_synth, chi_labels   = include_lobe(chi_synth, 0.2, 0.5)
    chi_synth_pv            = partial_volume_effect(chi_synth, 0.5, exclusions = exclusions, tissue_masks = chi_labels)
else
    chi_synth_pv            = partial_volume_effect(chi_synth, 0.5, exclusions = exclusions)
end

# Combine susceptibility sources
chi_local               = demean(chi_synth_pv, chi_mask)
chi_total               = chi_local + chi_background

# Compute fieldmaps (convolve with dipole kernel)
chi_local_padded        = anatomical_padding(chi_local)
chi_total_padded        = anatomical_padding(chi_total)
phi_local               = dipole_convolution(chi_local_padded)
phi_total               = dipole_convolution(chi_total_padded)
phi_local               = anatomical_unpadding(phi_local)
phi_total               = anatomical_unpadding(phi_total)

# ================ Save NIfTI files ================
save_nii(chi_mask, output_folder, "roi_mask.nii")              # ROI mask
save_nii(chi_local, output_folder, "chi_gt.nii")               # Local susceptibility (for dipole inversion evaluation)
save_nii(phi_local, output_folder, "phi_local_gt.nii")         # Local fieldmap (for background removal evaluation)
save_nii(phi_total, output_folder, "phi_total_gt.nii")         # Total fieldmap (for phase unwrapping evaluation)

println("Fieldmap simulation complete. Results saved to: ", output_folder)
# ===============================================