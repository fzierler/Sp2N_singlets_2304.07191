using HiRepAnalysis
using DelimitedFiles
using HDF5
# 
function Q_of_meas(conf_meas,conf,Qint)
    N = min(length(conf_meas),length(conf))
    Qs  = zeros(N)
    for ind in 1:N
        pos = findfirst(isequal(conf_meas[ind]),conf)
        Qs[ind] = isnothing(pos) ? NaN : Qint[pos]
    end
    return Qs
end
# read and parse both singlet and wilson_flow data
function Qs_and_index(fileQ,hdf5file,hdf5group)
    dat  = readdlm(fileQ,',',skipstart=1)
    conf = Int.(dat[:,1])
    Q    = dat[:,2]
    Qint = round.(Q)
    confh5 = HiRepAnalysis.h5property(hdf5file,hdf5group,"configurations")
    conf_meas = parse.(Int,last.(split.(confh5,'n')))

    # get all the indices where we have configurations of a given topological charge
    Qom = Q_of_meas(conf_meas,conf,Qint)
    voQ = filter(isfinite,unique(Qom))
    ind = [ findall(isequal(voQ[i]),Qom) for i in eachindex(voQ)]
    return voQ, ind
end
function splithdf5_Q(hdf5file,hdf5group,type,Qs,ind;absQ=false)
    hdf5f = h5open(hdf5file, "r")[hdf5group]
    measurements = filter( x -> contains(x,type), keys(hdf5f))
    for i in eachindex(ind)
        absQ && Qs[i] <= 0 && continue
        # create hdf file for fixed topological sector 
        # obtain path of the current
        if absQ
            file = "output/data_fixed_absQ.hdf5"
            if -Qs[i] ∈ Qs
                i_neg = findfirst(isequal(-Qs[i]),Qs)
                index = vcat(ind[i],ind[i_neg]) 
            else
                index = ind[i]
            end 
        else
            file = "output/data_fixed_Q.hdf5"
            index = ind[i]
        end
        hdf5groupQ = joinpath(hdf5group,"Q$((Qs[i]))")
        fid  = h5open(file,"cw")
        # loop over all unchanged entries measurements and merge
        identical = ["lofgile","sources","beta","gauge group","lattice","quarkmasses"]
        for property in identical
            if haskey(hdf5f,property)
                dat = read(hdf5f,property)
                write(fid,joinpath(hdf5groupQ,property),dat)
            end
        end
        for property in ["plaquette","configurations"]
            dat = read(hdf5f,property)
            write(fid,joinpath(hdf5groupQ,property),dat[index])
        end
        for Γ in measurements
            dat = read(hdf5f,Γ)
            # different dimensions for connected and disconnected measurements
            if ndims(dat) == 3
                write(fid,joinpath(hdf5groupQ,Γ),dat[index,:,:])
            elseif ndims(dat) == 2
                write(fid,joinpath(hdf5groupQ,Γ),dat[:,index])
            end
        end
        close(fid)
    end     
end
function mesons_in_flow_scale(hdf5file)
    Qdir = "output/Qhistories/"
    @showprogress for ensemble in readdir(Qdir)
        file = joinpath("output","Qhistories",ensemble)
        hdf5groupC = joinpath("runsSp4",ensemble,"out_spectrum")
        hdf5groupD = joinpath("runsSp4",ensemble,"out_spectrum_discon")
        
        typeD = "DISCON_SEMWALL SINGLET"
        typeC = "DEFAULT_SEMWALL TRIPLET"

        for absQ in [false,true]
            voQ,ind = Qs_and_index(file,hdf5file,hdf5groupD)
            splithdf5_Q(hdf5file,hdf5groupD,typeD,voQ,ind;absQ)
            voQ,ind = Qs_and_index(file,hdf5file,hdf5groupC)
            splithdf5_Q(hdf5file,hdf5groupC,typeC,voQ,ind;absQ)
        end
    end

    files = readdlm("input/parameters/param_deriv.csv",';',skipstart=1)[9:end,1]
    files = basename.(files)
    # flow scale w0 encoded in first line
    first_lines = readline.(joinpath.(Ref(Qdir),files))
    w0strings = getindex.(split.(first_lines,Ref(('(',')','='))),3)
    w0_with_uncertainty = split.(w0strings,Ref("+/-"))
    w0  = parse.(Float64,getindex.(w0_with_uncertainty,1))
    # combine this data with results from analysis in the julia scripts
    data,header = readdlm("output/data/data_deriv.csv",';',header=true)
    Δw0 = parse.(Float64,getindex.(w0_with_uncertainty,2))
    # remove all SU(2) simulations because we did not perform the gradient flow measurements there
    # we use that all SU(2) configurations have β=2.0
    data = data[9:end,:]
    # now combine results: 1st header, then data
    new_header = reshape(append!(vec(header),("w0","Deltaw0")),(1,29))
    new_data   = hcat(data,w0,Δw0)
    data_with_header = vcat(new_header,new_data)
    writedlm("output/data_flow.csv",data_with_header,';')
end
mesons_in_flow_scale(hdf5file)