using SparseArrays
using UnicodePlots

function read_csc_matrix(col_ptr_file::String, row_indices_file::String, values_file::String)
    # Read the column pointers, row indices, and values from files
    col_ptr = [parse(Int, line) for line in readlines(col_ptr_file)]
    row_indices = [parse(Int, line) for line in readlines(row_indices_file)]
    values = [parse(Float64, line) for line in readlines(values_file)]

    # Determine the number of rows and columns
    num_cols = length(col_ptr) - 1
    num_rows = maximum(row_indices) + 1  # Assuming 0-based indexing

    # Create the SparseMatrixCSC matrix
    A = SparseMatrixCSC{Float64, Int32}

    A->m = num_rows
    A->n = num_cols
    A->colptr = col_ptr
    A->rowval = row_indices
    A->nzval = values

    return A
end

let

  M = read_csc_matrix("mats/M_col_ptr", "mats/M_rows", "mats/M_vals")
  S = read_csc_matrix("mats/S_col_ptr", "mats/S_rows", "mats/S_vals")

  nothing

end
