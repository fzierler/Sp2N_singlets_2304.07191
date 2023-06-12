using HiRepAnalysis
using HDF5
using DelimitedFiles
function writehdf5_spectrum_disconnected(file,type,nhits;filterkey=false,key_pattern="",fnid="",only_channels=nothing,corronly=false,abspath="")
    filename = joinpath(abspath,joinpath(splitpath(file)[end-2:end]))
    path = joinpath(abspath,joinpath(splitpath(file)[end-2:end-1]))
    ispath(path) || mkpath(path)
    hdf5fn = filename*fnid*".hdf5"
    isfile(hdf5fn) && rm(hdf5fn)
    # save other relevant quantities
    if !corronly
        h5write(hdf5fn,"plaquette",plaquettes(file))
        h5write(hdf5fn,"configurations",confignames(file))
        h5write(hdf5fn,"gauge group",gaugegroup(file))
        h5write(hdf5fn,"beta",couplingβ(file))
        h5write(hdf5fn,"quarkmasses",quarkmasses(file))
        h5write(hdf5fn,"lattice",latticesize(file))
        h5write(hdf5fn,"lofgile",filename)
        h5write(hdf5fn,"sources",nhits)
    end
    # read correlator data
    c = correlators(file,type;withsource=true,average=false,filterkey,key_pattern)
    if isnothing(only_channels)
        channels = unique(collect(keys(c[1])))
    else
        channels = intersect(unique(collect(keys(c[1]))),only_channels)
    end
    for Γ in channels
        d = getindex.(c,Γ)
        dat = flatten_disc(d,nhits;rescale=1)
        h5write(hdf5fn,type*"_"*Γ,dat)
    end
end
function writehdf5_spectrum_connected(file,type;abspath="",fnid = "")
    filename = joinpath(abspath,joinpath(splitpath(file)[end-2:end]))
    path = joinpath(abspath,joinpath(splitpath(file)[end-2:end-1]))
    ispath(path) || mkpath(path)
    hdf5fn = filename*fnid*".hdf5"
    isfile(hdf5fn) && rm(hdf5fn)
    # save other relevant quantities
    h5write(hdf5fn,"plaquette",plaquettes(file))
    h5write(hdf5fn,"configurations",confignames(file))
    h5write(hdf5fn,"gauge group",gaugegroup(file))
    h5write(hdf5fn,"beta",couplingβ(file))
    h5write(hdf5fn,"quarkmasses",quarkmasses(file))
    h5write(hdf5fn,"lattice",latticesize(file))
    h5write(hdf5fn,"lofgile",filename)
    h5write(hdf5fn,"measurement_type",type)
    # read correlator data
    c = correlators(file,type;withsource=false,average=true)
    channels = unique(collect(keys(c[1])))
    for Γ in channels
        dat = reduce(hcat,getindex.(c,Γ))
        h5write(hdf5fn,type*"_"*Γ,dat)
    end
end
function writehdf5_spectrum_connected(file,typeU,typeD,typeUD;abspath="")
    fnid = ""
    # filename
    filename = joinpath(abspath,joinpath(splitpath(file)[end-2:end]))
    path = joinpath(abspath,joinpath(splitpath(file)[end-2:end-1]))
    ispath(path) || mkpath(path)
    hdf5fn = filename*fnid*".hdf5"
    isfile(hdf5fn) && rm(hdf5fn)
    # save other relevant quantities
    h5write(hdf5fn,"plaquette",plaquettes(file))
    h5write(hdf5fn,"configurations",confignames(file))
    h5write(hdf5fn,"gauge group",gaugegroup(file))
    h5write(hdf5fn,"beta",couplingβ(file))
    h5write(hdf5fn,"quarkmasses",quarkmasses(file))
    h5write(hdf5fn,"lattice",latticesize(file))
    h5write(hdf5fn,"lofgile",filename)
    h5write(hdf5fn,"measurement_type",[typeU,typeD,typeUD])
    # read correlator data
    for type in (typeU,typeD,typeUD)
        c = correlators(file,type;withsource=false,average=true)
        channels = unique(collect(keys(c[1])))
        for Γ in channels
            dat = getindex.(c,Γ)
            h5write(hdf5fn,type*"_"*Γ,reduce(hcat,dat))
        end
    end
end
function writehdf5_spectrum_disconnected_nondeg(file,type,nhits,masses;filterkey=false,key_pattern="",fnid="",corronly=false,abspath="",only_channels=nothing)
    # filename
    filename = joinpath(abspath,joinpath(splitpath(file)[end-2:end]))
    path = joinpath(abspath,joinpath(splitpath(file)[end-2:end-1]))
    ispath(path) || mkpath(path)
    hdf5fn = filename*fnid*".hdf5"
    isfile(hdf5fn) && rm(hdf5fn)
    # save other relevant quantities
    if !corronly
        h5write(hdf5fn,"plaquette",plaquettes(file))
        h5write(hdf5fn,"configurations",confignames(file))
        h5write(hdf5fn,"gauge group",gaugegroup(file))
        h5write(hdf5fn,"beta",couplingβ(file))
        h5write(hdf5fn,"quarkmasses",quarkmasses(file))
        h5write(hdf5fn,"lattice",latticesize(file))
        h5write(hdf5fn,"lofgile",filename)
        h5write(hdf5fn,"sources",nhits)
    end
    # read correlator data
    c1 = correlators(file,type;withsource=true,average=false,masses=true,mass=masses[1],filterkey,key_pattern)
    c2 = correlators(file,type;withsource=true,average=false,masses=true,mass=masses[2],filterkey,key_pattern)
    if isnothing(only_channels)
        channels = unique(collect(keys(c1[1])))
    else
        channels = intersect(unique(collect(keys(c1[1]))),only_channels)
    end
    for Γ in channels
        d1 = getindex.(c1,Γ)
        d2 = getindex.(c2,Γ)
        dat1 = flatten_disc(d1,nhits;rescale=1)
        dat2 = flatten_disc(d2,nhits;rescale=1)
        h5write(hdf5fn,type*"_"*Γ*"_"*masses[1],dat1)
        h5write(hdf5fn,type*"_"*Γ*"_"*masses[2],dat2)
    end
end
function mergehdf_discon(file,type,file1,file2)
    isfile(file) && rm(file)
    hdf5_f1 = h5open(file1, "r")
    hdf5_f2 = h5open(file2, "r")
    # test that we used the same sources in both cases
    @assert keys(hdf5_f1) == keys(hdf5_f2)     
    disc_measurements = filter( x -> contains(x,type), keys(hdf5_f1))
    # test that we compare data from the same lattice setup
    identical = ["beta","gauge group", "lattice",  "plaquette" ,"quarkmasses","configurations"]
    for property in identical
        @assert read(hdf5_f1,property) == read(hdf5_f2,property) "$property mismatch"
        h5write(file,property,read(hdf5_f1,property))
    end
    # loop over all measurements and merge
    N1 = size(read(hdf5_f1,disc_measurements[1]))[1]
    N2 = size(read(hdf5_f2,disc_measurements[1]))[1]
    if N1 != N2 
        @warn "mismatch file1 N=$N1 and file N=$N2"
    end
    N = min(N1,N2)
    for Γ in disc_measurements
        dat = cat(read(hdf5_f1,Γ)[1:N,:,:],read(hdf5_f2,Γ)[1:N,:,:],dims=2)
        h5write(file,Γ,dat)
    end
    h5write(file,"lofgile",file1*file1)
    h5write(file,"sources",read(hdf5_f1,"sources")+read(hdf5_f2,"sources"))
end
function mergehdf_discon(file,type,file1)
    isfile(file) && rm(file)
    cp(file1,file;force=true)
end
function fermion_masses_from_filename(file)
    p1 = last(findfirst("m1",file)) 
    p2 = first(findfirst("m2",file))
    p3 = last(findfirst("m2",file)) 
    p4 = findnext('/',file,p3) 
    # create array of masses for matching of output
    m  = [file[p1+1:p2-1],file[p3+1:p4-1]] 
    return m
end
function logfiles_to_hdf5(name,path,hdfpath,fileC,typeC,typeD)
    dir = joinpath(path,name) 
    discon_logs = filter(contains("discon"),readdir(dir,join=true))
    # the fifth element after splitting contains the number of hits
    hits_strings = getindex.(split.(basename.(discon_logs),"_"),4)
    discon_hits  = parse.(Int,replace.(hits_strings,"h"=>""))
    # restrict ourselves to scalar and pseudoscalar data
    only_channels = ["id_disc_re", "g5_disc_re"]
    writehdf5_spectrum_connected(fileC,typeC;abspath=hdfpath)
    for (log,hits) in zip(discon_logs,discon_hits)
        writehdf5_spectrum_disconnected(log,typeD,hits;abspath=hdfpath,only_channels)
    end
end
# merge all generated hdf5 files for disconnected contributions
function merge_all_hdf5_files(name,hdfpath,typeD)
    dir = joinpath(hdfpath,name)
    discon_hdf5 = filter(contains("discon"),readdir(dir,join=true))
    tmpname = joinpath(dir,"tmp.hdf5")
    dstname = joinpath(dir,"out_spectrum_discon.hdf5")
    for (i,log) in enumerate(discon_hdf5)
        if i==1
            if length(discon_hdf5) > 1
                mergehdf_discon(tmpname,typeD,log)
            else
                mergehdf_discon(dstname,typeD,log)
            end
        else
            mergehdf_discon(dstname,typeD,log,tmpname)
            cp(dstname,tmpname;force=true)
        end
    end
    isfile(tmpname) && rm(tmpname)
    rm.(discon_hdf5)
end
function single_hdf5(hdfpath,new_file)
    hdf5 = h5open(new_file,"w")
    for runs in readdir(hdfpath;join=true)
        for ensemble in readdir(runs;join=true)
            for file in readdir(ensemble;join=true)
                h5file = relpath(file,joinpath(hdfpath))
                h5base = splitext(h5file)[1]
                f = h5open(file,"r")
                for k in keys(f)
                    handle = joinpath(h5base,k)
                    data = read(f,k)
                    write(hdf5,handle,data)
                end
            end
        end
    end
end
function write_raw_to_hdf5(path,hdfpath,hdf5file)
    namesDeg    = readdlm("input/parameters/param_both.csv",';';skipstart=1)[:,1]
    namesNonDeg = readdlm("input/parameters/param_non_deg.csv",';';skipstart=1)[:,1]
    hitsNonDeg  = readdlm("input/parameters/param_non_deg.csv",';';skipstart=1)[:,2]
    type   = "DISCON_SEMWALL SINGLET"
    typeC  = "DEFAULT_SEMWALL TRIPLET"
    typeU  = "SEMWALL_U TRIPLET"
    typeD  = "SEMWALL_D TRIPLET"
    typeUD = "SEMWALL_UD TRIPLET"
    filesC1 = joinpath.(path,namesDeg,"out_spectrum")
    filesC2 = joinpath.(path,namesNonDeg,"out_spectrum_with_pcac")

    # special case degenerate limit for non-degenerate data
    writehdf5_spectrum_connected(filesC2[8],typeC;abspath=hdfpath)
    # special case smeared analysis output file
    smeared_file = joinpath(path,"runsSp4/Lt32Ls16beta6.9m1-0.90m2-0.90/out_spectrum_smeared")
    typeS1 = "source_N60_sink_N0 TRIPLET"
    typeS2 = "source_N60_sink_N60 TRIPLET"
    writehdf5_spectrum_connected(smeared_file,typeS1;abspath=hdfpath,fnid="_N60N0")
    writehdf5_spectrum_connected(smeared_file,typeS2;abspath=hdfpath,fnid="_N60N60")

    # parse non-degenerate raw logfiles
    println("raw log files with non-degenerate fermion masses...")
    @showprogress for i in eachindex(namesNonDeg[1:7])
        hits  = hitsNonDeg[i]
        fileD = joinpath(path,namesNonDeg[i],"out_spectrum_discon_h$hits")
        fileDtmp = joinpath(hdfpath,namesNonDeg[i],"out_spectrum_discon_h$hits.hdf5")
        fileDout = joinpath(hdfpath,namesNonDeg[i],"out_spectrum_discon.hdf5")
        masses = fermion_masses_from_filename(fileD)
        writehdf5_spectrum_disconnected_nondeg(fileD,type,hits,masses;abspath=hdfpath,only_channels=["id_disc_re", "g5_disc_re"])
        mv(fileDtmp,fileDout)
        writehdf5_spectrum_connected(filesC2[i],typeU,typeD,typeUD;abspath=hdfpath)
    end
    # parse degenerate raw logfiles
    println("raw log files with degenerate fermion masses...")
    @showprogress for (i,name) in enumerate(namesDeg)
        logfiles_to_hdf5(name,path,hdfpath,filesC1[i],typeC,type)
        merge_all_hdf5_files(name,hdfpath,type)
    end
    # Combine individual files into a single hdf5 file
    single_hdf5(hdfpath,hdf5file)
end
write_raw_to_hdf5(path,hdfpath,hdf5file)