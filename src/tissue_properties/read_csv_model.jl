"""
    read_csv_model(model_file::String) -> Matrix{Float64}

Reads a `.csv` file containing tissue model parameters and converts 
it into a numerical matrix.

# Arguments
- `model_file::String`: Path to the `.csv` file containing the model. The file is 
  assumed to have one row per tissue and at least one column of metadata 
  followed by numerical parameters.

# Returns
- `model::Matrix{Float64}`: A numeric matrix of size (N × M), where `N` is the 
  number of tissues and `M` is the number of model parameters. The first column 
  of the CSV is assumed to be metadata and is skipped.

# Example
    chi_model = read_csv_model("path/to/chi_model.csv")

# Notes
- Only columns 2 to end are considered numeric and extracted.
- The CSV must be formatted such that all non-numeric fields are in the first column.


Last edited by Silva J on 2025.08.03. 
"""
function read_csv_model(model_file::String)    
    # Read csv files
    model_df    = CSV.read(model_file,DataFrame)

    # Convert into matrix with numerical values
    model       = Matrix(model_df[:, 2:end])

    return model
end