# =========================================================================================
# example01_synthetic_chi.jl
# =========================================================================================
#
# Example script for generating synthetic susceptibility maps using QSMPhantom.
#
# To run this script from the project root directory:
# - In Julia REPL: 
#       include("examples/example01_synthetic_chi.jl")
# - From the command line: 
#       julia --project=. -e 'include("examples/example01_synthetic_chi.jl")'
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
using QSMPhantom

# ========= User-adjustable parameters ==========
output_folder       = joinpath(@__DIR__, "..", "results", "ex01_synthetic_chi")
case_select         = 1             # 1 = healthy_subject; 2 = iron_overload; 3 = pathological_lobe
exclusions          = [12,21,23]    # Exclude heart (12), large intestine (21), and lungs (23)

# =================== MAIN CODE =================
# Select susceptibility model
cases                   = ["healthy_subject", "iron_overload","pathological_lobe"]
chi_model_file          = joinpath(@__DIR__, "..", "src", "models", cases[case_select], "chi_model.csv")
chi_model               = read_csv_model(chi_model_file)

# Compute local susceptibility
chi_synth, chi_mask     = local_susceptibility_per_tissue(chi_model)

#Compute background susceptibilities
air_chi_ppm             = 9.2   # Susceptibility value for external air (in ppm)
chi_synth, chi_mask     = local_susceptibility_exclude(chi_synth, chi_mask, exclusions)
chi_background          = background_susceptibility(chi_model, air_chi_ppm, exclusions)

# Add lobe (if selected) and apply partial volume effects
if case_select == 3
    chi_synth, chi_labels   = include_lobe(chi_synth, 0.2, 0.5)
    chi_synth_pv            = partial_volume_effect(chi_synth, 0.5, exclusions = exclusions, tissue_masks = chi_labels)
else
    chi_synth_pv            = partial_volume_effect(chi_synth, 0.5, exclusions = exclusions)
end

#Demean local susceptibility and compute total susceptibility
chi_local               = demean(chi_synth_pv,chi_mask)
chi_total               = chi_local + chi_background

# Save output files
save_nii(chi_mask, output_folder, "roi_mask.nii")               # ROI mask  
save_nii(chi_local, output_folder, "chi_local.nii")             # Local susceptibility
save_nii(chi_background, output_folder, "chi_background.nii")   # Background susceptibility
save_nii(chi_total, output_folder, "chi_total.nii")             # Total susceptibility

println("Simulation complete. Results saved to: ", output_folder)
# ===============================================