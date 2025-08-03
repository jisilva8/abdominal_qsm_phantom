# =========================================================================================
# example04_simulate_acquisition.jl
# =========================================================================================
#
# Example script for simulating complex water+fat signal acquisitions using QSMPhantom.
#
# This simulation includes ground truth maps for susceptibility, fieldmaps, R2*, PDFF,
# and corresponding complex water/fat signal components. The generated dataset is 
# suitable for validating water-fat separation, background field removal, and QSM methods.
#
# To run this script from the project root directory:
# - In Julia REPL:
#       include("examples/example04_simulate_acquisition.jl")
# - From the command line:
#       julia --project=. -e 'include("examples/example04_simulate_acquisition.jl")'
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
output_folder       = joinpath(@__DIR__, "..", "results", "ex04_simulate_acquisition")
case_select         = 1             # 1 = healthy_subject; 2 = iron_overload; 3 = pathological_lobe
exclusions          = [12, 21, 23]  # Exclude heart (12), large intestine (21), and lungs (23)
Npeaks              = 6             # Fat model: number of peaks (options: 1,3,4,5,6,7,9)
rBW                 = 400.0         # Receiver bandwidth in Hz/px (e.g., 400 or 2000)
TE_select           = 1             # Echo time configuration:
                                    # 1 = 5-echo in-phase echo times (single-peak)
                                    # 2 = 5-echo effective in-phase echo times (multi-peak)
                                    # 3 = 6-echo CSE acquisition for graph cut methods

# =================== MAIN CODE =================
# === Define cases and echo times ===
cases                   = ["healthy_subject", "iron_overload", "pathological_lobe"]
TE_inphase              = [2.30, 4.60, 6.90, 9.20, 11.50]            # Single-peak in-phase TE
TE_effinphase           = [2.34, 4.57, 6.77, 9.14, 11.59]            # Multi-peak effective in-phase TE
TE_cse                  = [1.8, 3.36, 4.92, 6.48, 8.04, 9.6]         # CSE acquisition
TE_ms                   = [TE_inphase, TE_effinphase, TE_cse]

# === Simulate Susceptibility map ===
# Select susceptibility model
chi_model_file          = joinpath(@__DIR__, "..", "src", "models", cases[case_select], "chi_model.csv")
chi_model               = read_csv_model(chi_model_file)

# Compute local susceptibility
chi_synth, chi_mask     = local_susceptibility_per_tissue(chi_model)

# Compute background susceptibilities
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

# Demean local susceptibility and compute total susceptibility
chi_local               = demean(chi_synth_pv,chi_mask)
chi_total               = chi_local + chi_background

# === Simulate Fieldmap ===
# Compute fieldmaps (convolve with dipole kernel)
chi_local_padded        = anatomical_padding(chi_local)
chi_total_padded        = anatomical_padding(chi_total)
phi_local               = dipole_convolution(chi_local_padded)
phi_total               = dipole_convolution(chi_total_padded)
phi_local               = anatomical_unpadding(phi_local)
phi_total               = anatomical_unpadding(phi_total)


# === Simulate R2* map ===
# Select R2* model file
r2star_model_file           = joinpath(@__DIR__, "..", "src", "models", cases[case_select], "r2s_model.csv")

# Generate synthetic R2* map
r2star_model                = read_csv_model(r2star_model_file)
r2star_synth, r2star_labels = r2star_per_tissue(r2star_model)

# Optionally include lobe with contrast
if case_select == 3
    r2star_synth, r2star_labels = include_lobe(r2star_synth, 150.0, 850.0, tissue_masks = r2star_labels)
end

# Apply partial volume effects
r2star_synth                = partial_volume_effect(r2star_synth, 0.8, tissue_masks = r2star_labels)


# === Acquisition Sequence ===
# Simulate multiecho signal
signal, water, fat, pdff    = simulate_acquisition(phi_total, r2star_synth, TE_ms[TE_select], chi_mask, Npeaks, rBW)


# === Save output files ===
save_nii(chi_mask, output_folder, "roi_mask.nii")           # ROI mask
save_nii(chi_local, output_folder, "chi_gt.nii")            # Local susceptibility ground truth
save_nii(phi_local, output_folder, "phi_local_gt.nii")      # Local fieldmap (for background removal)
save_nii(phi_total, output_folder, "phi_total_gt.nii")      # Total fieldmap (for phase unwrapping)
save_nii(r2star_synth, output_folder, "r2star.nii")         # R2* map
save_nii(pdff, output_folder, "pdff_gt.nii")                # Proton density fat fraction ground truth

save_nii(abs.(signal), output_folder, "signal_mag.nii")     # Total signal magnitude (water + fat)
save_nii(angle.(signal), output_folder, "signal_phs.nii")   # Total signal phase
save_nii(abs.(water), output_folder, "water_mag.nii")       # Water signal magnitude
save_nii(angle.(water), output_folder, "water_phs.nii")     # Water signal phase
save_nii(abs.(fat), output_folder, "fat_mag.nii")           # Fat signal magnitude
save_nii(angle.(fat), output_folder, "fat_phs.nii")         # Fat signal phase

println("Acquisition simulation complete. Results saved to: ", output_folder)
# ===============================================