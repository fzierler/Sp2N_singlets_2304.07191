# Perform analysis based on the hdf5 file
# and the directory containing the gradient flow data
hdf5file = "input_data/data.hdf5"
flowdir  = "input_data/wilson_flow/"

# Allows to generate the hdf5file from the raw logs
# path    specifies the directory containing the raw log files
# hdfpath specifies temporary directory for saving individual hdf5 files foir every ensemble
start_from_raw_logs = false
path_raw_logs = "../Singlets_Data/"
hdfpath  = "input_data/hdf_files/"
include("singlet_anaysis.jl")
