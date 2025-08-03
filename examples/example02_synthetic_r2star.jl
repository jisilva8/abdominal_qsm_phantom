# =========================================================================================
# example02_synthetic_r2star.jl
# =========================================================================================
#
# Example script for generating synthetic R2* (1/T2*) maps using QSMPhantom.
#
# To run this script from the project root directory:
# - In Julia REPL: 
#       include("examples/example02_synthetic_r2star.jl")
# - From the command line: 
#       julia --project=. -e 'include("examples/example02_synthetic_r2star.jl")'
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
output_folder       = joinpath(@__DIR__, "..", "results", "ex02_synthetic_r2star")
case_select         = 1             # 1 = healthy_subject; 2 = iron_overload; 3 = pathological_lobe

# =================== MAIN CODE =================
# Select R2* model file
cases                   = ["healthy_subject", "iron_overload", "pathological_lobe"]
r2star_model_file       = joinpath(@__DIR__, "..", "src", "models", cases[case_select], "r2s_model.csv")
r2star_model            = read_csv_model(r2star_model_file)

# Compute local susceptibility
r2star_synth, r2_labels = r2star_per_tissue(r2star_model)

# Optionally include lobe with contrast
if case_select == 3
    r2star_synth, r2_labels = include_lobe(r2star_synth, 150.0, 850.0, tissue_masks = r2_labels)
end

# Apply partial volume effects
r2star_synth = partial_volume_effect(r2star_synth, 0.8, tissue_masks = r2_labels)

# Save output file
save_nii(r2star_synth, output_folder, "r2star.nii")  # Synthetic R2* map

println("R2* simulation complete. Results saved to: ", output_folder)
# ===============================================