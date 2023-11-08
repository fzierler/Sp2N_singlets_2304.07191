using HDF5
using DelimitedFiles
using LsqFit
using Statistics
using LinearAlgebra
include("variational_analysis.jl")

h5dir = "./output/matrix_correlators/"
hdf5files = readdir(h5dir,join=true)

for file in hdf5files
    # set up directory for saving eigenvalues
    eigval_file = joinpath("output","eigvals",basename(file)) 
    ispath(dirname(eigval_file)) || mkpath(dirname(eigval_file)) 
   
    # obtain eigenvalues from a jackknife analysis
    file_id = h5open(file,"r")
   
    # solve GEVP
    corr1 = read(file_id,"correlator_matrix_A") 
    corr1_deriv = read(file_id,"correlator_matrix_A_deriv") 
    eigvals, Δeigvals = eigenvalues_jackknife(corr1)
    eigvals_deriv, Δeigvals_deriv = eigenvalues_jackknife(corr1_deriv)

    # create a new h5 file
    file_id_eigval = h5open(eigval_file,"w")
    # save data that has been unchanged
    for k in ["beta","configurations","gauge group","lattice","logfile","plaquette","quarkmasses"]
        file_id_eigval[k] = read(file_id,k) 
    end

    # write eigenvalues
    file_id_eigval["eigvals"] = eigvals
    file_id_eigval["Δeigvals"] = Δeigvals
    file_id_eigval["eigvals_deriv"] = eigvals_deriv
    file_id_eigval["Δeigvals_deriv"] = Δeigvals_deriv

    # close files 
    close(file_id)
    close(file_id_eigval)
end