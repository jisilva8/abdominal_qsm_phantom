"""
    save_nii(Im::Union{Array{Float64,3}, Array{Float64,4}}, folder::String, filename::String) -> Nothing

Saves a 3D or 4D array as a NIfTI file in the specified folder with the given filename.

# Arguments
- `Im::Union{Array{Float64,3}, Array{Float64,4}}`: 3D or 4D array representing the image volume(s).
- `folder::String`: Path to the folder where the NIfTI file will be saved. The folder is created if it does not exist.
- `filename::String`: Name of the output NIfTI file.

# Returns
- `Nothing`: This function saves the file to disk and prints a confirmation message.

# Example
    save_nii(image_volume, "output/nifti", "subject01.nii")

# Notes
- The function checks that the input array is 3D or 4D.
- Uses `NIVolume` and `niwrite` for conversion and saving.

Last edited by Silva J on 2025.08.03.
"""
function save_nii(Im::Union{Array{Float64,3}, Array{Float64,4}}, folder::String, filename::String)
    #Check dimensionality
    ndims(Im) in (3, 4) || error("Input array must be 3D or 4D. Got $(ndims(Im))D.")

    #Create folder if it does not exist
    isdir(folder) || mkpath(folder)

    # Full and relative file path
    filepath = joinpath(folder, filename)
    relative_path = relpath(filepath, pwd())

    #Convert to NIfTI and save
    volume = NIVolume(Im)
    niwrite(filepath, volume)

    println("✅ Saved NIfTI to: ", relative_path)
end