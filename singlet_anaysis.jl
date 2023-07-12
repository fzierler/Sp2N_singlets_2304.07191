# Install required packages
using Pkg
Pkg.activate(".")
Pkg.instantiate()

# include required Packages
using ProgressMeter

path = path_raw_logs
isdir("output") || mkdir("output")
# run python scripts for performing analysis of
println("Gradient flow analysis...")
# this is not very elegant
# modify python script to modify the path of the wilso_flow directory
function change_flowdir(flowdir,py_script;tmp="tmp_flow_analysis.py")
    fid = open(tmp,"w+")
    modified = false
    for line in eachline(py_script)
        if startswith(line,"path") && !modified
            write(fid,"path = '$flowdir'\n")
            modified = true
        else
            write(fid,line*"\n")
        end
    end
    close(fid)
    mv(tmp,py_script;force=true)
end
change_flowdir(flowdir,"flow_analysis.py")
# The script assumes that python or python3 is available and
# the rquired packages are installed
try
    run(`python flow_analysis.py`)
catch
    run(`python3 flow_analysis.py`)
end

if start_from_raw_logs
    # first convert all log files to hdf5
    isdir(hdfpath) || mkdir(hdfpath)
    include("scripts/writeHDF5.jl")
end

# perform analysis of degenerate and non-degenerate mesons
println("Analysis of degenerate fermion ensembles...")
include("scripts/singlet_deg.jl")
println("Analysis of non-degenerate fermion ensembles...")
include("scripts/singlet_nondeg.jl")
# then combine these results with the wilson flow results
# additionally split hdf5 files into files with fixed topological charge Q
println("Analysis of fixed topological sectors...")
include("scripts/Qhistory.jl")
# write tables
println("Write tables...")
include("scripts/singlet_tables.jl")
# create Plots
println("Generate figures...")
isdir("output/figures") || mkdir("output/figures")
include("scripts/singlet_plots/plot_comparison.jl")
include("scripts/singlet_plots/plot_singlet_constant.jl")
include("scripts/singlet_plots/plot_singlet_flow.jl")
include("scripts/singlet_plots/plot_singlet_nondeg.jl")
include("scripts/singlet_plots/plot_singlet.jl")

# perform analysis with smeared correlators
println("Analysis of smeared operator...")
include("scripts/analysis_smeared.jl")
println("done.")
