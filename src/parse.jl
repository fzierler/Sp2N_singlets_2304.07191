function gaugegroup(file)
    if HDF5.ishdf5(file)
        hdf5 = h5open(file, "r")
        return read(hdf5,"gauge group")
    else
        return gaugegroup_log(file)
    end
end
function quarkmasses(file)
    if HDF5.ishdf5(file)
        hdf5 = h5open(file, "r")
        return read(hdf5,"quarkmasses")
    else
        return quarkmasses_log(file)
    end
end
function latticesize(file)
    if HDF5.ishdf5(file)
        hdf5 = h5open(file, "r")
        return read(hdf5,"lattice")
    else
        return latticesize_log(file)
    end
end
function correlators(file,type,key;kws...)
    if HDF5.ishdf5(file)
        hdf5 = h5open(file, "r")
        return read(hdf5,"$(type)_$(key)")
    else
        return correlators_logfile(file,type,key;kws...)
    end
end
function correlators(file,type,key,nhits::Int;maxhits,masses=false,mass=nothing,kws...)
    if HDF5.ishdf5(file)
        hdf5 = h5open(file, "r")
        if !masses
           corr = read(hdf5,"$(type)_$(key)")[:,1:min(maxhits,nhits),:]
        else
           corr = read(hdf5,"$(type)_$(key)_$(mass)")[:,1:min(maxhits,nhits),:]
        end
        return corr
    else
        d = correlators_logfile(file,type,key;masses,mass,kws...)
        return flatten_disc(d,nhits;rescale=1)
    end
end
function plaquettes(file;therm=0,step=1)
    if HDF5.ishdf5(file)
        hdf5 = h5open(file, "r")
        return read(hdf5,"plaquette")[therm+1:step:end]
    else
        return plaquettes_log(file;therm,step)
    end
end
function nconfigs(file)
    if HDF5.ishdf5(file)
        hdf5 = h5open(file, "r")
        return length(read(hdf5,"configurations"))
    else
        return nconfigs_log(file)
    end
end
function gaugegroup_log(file)
    for line in eachline(file)
        if occursin("Gauge group",line)
            pos = findlast(' ',line)
            return strip(line[pos:end])
        end
    end
end
function latticesize_log(file)
    for line in eachline(file)
        if occursin("Global size is",line)
            pos  = last(findfirst("Global size is",line))+1
            sizestring  = lstrip(line[pos:end])
            latticesize = parse.(Int,split(sizestring,"x"))
            return latticesize
        end
    end
end
function correlators_logfile(file,type,key;kws...)
    corrs = correlators(file,type;filterkey=true,key_pattern=key,kws...)
    return reduce(hcat,getindex.(corrs,key))
end
function correlators(file,type;withsource=false,average=true,masses=false,mass="",filterkey=false,key_pattern="")
    T = latticesize(file)[1]
    corr = zeros(T) # preallocate array for parsing of correlator
    dict = Dict{String,Vector{Float64}}()
    dictarray = Dict{String,Vector{Float64}}[]
    conf0 = 0
    src0  = 0
    # keep track of position in file for progress meter
    for line in eachline(file)
        if occursin(type,line)
            if masses
                occursin("mass=$mass",line) || continue
            end
            if filterkey
                occursin(key_pattern,line) || continue
            end
            # get configuration number
            pos_num = findfirst('#',line)
            end_num = findnext(' ',line,pos_num)
            conf = parse(Int,line[pos_num+1:end_num-1])
            # find number of the source if available
            if withsource
                pos_src = last(findfirst("src",line))+1
                end_src = findnext(' ',line,pos_src+1)
                src = parse(Int,line[pos_src:end_src])
            else
                src = 0
            end
            # find last '=' sign which separates values from Γ structure
            # TODO this does not work for momenta
            pos_eq = findlast('=',line)
            #key_st = findprev(' ',line,pos_eq)
            key_st = last(findfirst(type,line))+1
            key = line[key_st+1:pos_eq-1]
            if withsource
                # create new entry if configuration or source number changes
                # if we need to parse more than one source at a time per configuration
                if conf0 != conf || src0 != src
                    if !isempty(dict)
                        push!(dictarray,dict)
                        dict = Dict{String,Vector{Float64}}()
                    end
                end
            end
            # parse corrrelator values
            pos_0 = findnext(' ',line,pos_eq)
            for t in 1:T
                pos_1 = findnext(' ',line,pos_0+1)
                corr[t] = Parsers.parse(Float64,line[pos_0:pos_1])
                pos_0 = pos_1
            end
            dict[key] = copy(corr)
            conf0 = conf
            src0  = src
        end
        if !withsource
            # If we only have one source at a time and possibly one configuration
            # at a time. Under these circumstances the method used to separate distinct
            # measurements fails. In this case the end of measurement on a given confiuration
            # is signalled by a line that reads:
            # [MAIN][0]Configuration #N: analysed in [a sec b usec]
            if occursin("analysed",line)
                if !isempty(dict)
                    push!(dictarray,dict)
                    dict = Dict{String,Vector{Float64}}()
                end
            end
        end
    end
    if !isempty(dict)
        push!(dictarray,dict)
    end
    if average
        averagevectors!(dictarray,T)
    end
    return dictarray
end
function plaquettes_log(file;therm=0,step=1)
    plaquettes = Float64[]
    for line in eachline(file)
        if occursin("Plaquette",line)
            line = replace(line,"="=>" ")
            line = replace(line,":"=>" ")
            p = parse(Float64,split(line)[end])
            append!(plaquettes,p)
        end
    end
    return plaquettes[therm+1:step:end]
end
function HMC_accept(file)
    acceptreject = Bool[]
    ΔS = Float64[]
    for line in eachline(file)
        if occursin("DeltaS",line)
            pos1 = last(findfirst("DeltaS",line))+1
            pos2 = first(findfirst("exp(-DS)",line))-1
            Sstring = line[pos1:pos2]
            for symb in ["=","[","]"]
                Sstring = replace(Sstring,symb=>" ")
            end
            append!(ΔS,parse(Float64,Sstring))
        end
        if occursin("Configuration rejected.",line)
            append!(acceptreject,false)
        end
        if occursin("Configuration accepted.",line)
            append!(acceptreject,true)
        end
    end
    return BitArray(acceptreject), ΔS
end
function loweig(file,therm=0)
    eigs = Float64[]
    for line in eachline(file)
        if occursin("[LOWEIG][0]Eig 0 =",line)
            eig = parse(Float64,split(line,"=")[end])
            append!(eigs,eig)
        end
    end
    return eigs[therm+1:end]
end
function quarkmasses_log(file)
    masses = Float64[]
    for line in eachline(file)
        if occursin("[MAIN][0]Mass[0]",line)
            s = split(line,(","))
            for i in eachindex(s)
                m = parse(Float64,split(s[i],"=")[2])
                append!(masses,m)
            end
            return masses
        end
    end
end
function couplingβ(file)
    try
        l = split(file,"beta")[end]
        β = parse(Float64,split(l,"m")[1])
        return β
    catch
        for line in eachline(file)
            if occursin("Configuration from",line)
                l = split(line,"b")[end]
                l = split(l,"m")[1]
                β = parse(Float64,l)
                return β
            end
        end
    end
end
function fileinfo(file)
    Λ = latticesize(file)
    q = quarkmasses(file)
    β = couplingβ(file)
    return Λ,q,β
end
#################################################
# Disconnected Measurements from /Disocnnected  #
#################################################
function dilution(file)
    for line in eachline(file)
        if occursin("will be used",line)
            eo    = occursin("eo"   ,lowercase(line))
            time  = occursin("time" ,lowercase(line))
            spin  = occursin("spin" ,lowercase(line))
            color = occursin("color",lowercase(line))
            return eo, time, spin, color
        end
    end
end
function ncolors(file)
    for line in eachline(file)
        if occursin("Gauge group",line)
            pos1 = findfirst('(',line)+1
            pos2 = findfirst(')',line)-1
            colors = parse(Int,line[pos1:pos2])
            return colors
        end
    end
end
function nconfigs_log(file)
    nconfig = 0
    for line in eachline(file)
        if occursin("read",line)
            if occursin("Configuration",line)
                nconfig += 1
            end
        end
    end
    return nconfig
end
function nhits(file)
    for line in eachline(file)
        if occursin("Number of noise vector : nhits",line)
            pos = findfirst('=',line)
            nhits =  parse(Int,line[pos+1:end])
            return nhits
        end
    end
end
function confignames(file)
    fns = AbstractString[]
    for line in eachline(file)
        if occursin("read",line)
            if occursin("Configuration",line)
                pos1 = findlast('/',line)
                pos2 = findnext(']',line,pos1)
                push!(fns,line[pos1+1:pos2-1])
            end
        end
    end
    return fns
end
function h5property(hdf5file,hdf5group,property)
    @assert HDF5.ishdf5(hdf5file)
    hdf5 = h5open(hdf5file, "r")
    prop = joinpath(hdf5group,property)
    return read(hdf5,prop)
end
gaugegroup(hdf5file,hdf5group) = h5property(hdf5file,hdf5group,"gauge group")
latticesize(hdf5file,hdf5group) = h5property(hdf5file,hdf5group,"lattice")
quarkmasses(hdf5file,hdf5group) = h5property(hdf5file,hdf5group,"quarkmasses")
nconfigs(hdf5file,hdf5group) = length(h5property(hdf5file,hdf5group,"configurations"))
couplingβ(hdf5file,hdf5group) = h5property(hdf5file,hdf5group,"beta")
correlators(hdf5file,hdf5group,type,key;kws...) = h5property(hdf5file,hdf5group,"$(type)_$(key)")
plaquettes(hdf5file,hdf5group;therm=0,step=1) = h5property(hdf5file,hdf5group,"plaquette")[therm+1:step:end]
function correlators(hdf5file,hdf5group,type,key,nhits::Int;maxhits,masses=false,mass=nothing,kws...)
    @assert HDF5.ishdf5(hdf5file)
    hdf5 = h5open(hdf5file, "r")
    if !masses
        prop = joinpath(hdf5group,"$(type)_$(key)")
    else
        prop = joinpath(hdf5group,"$(type)_$(key)_$(mass)")
    end
    corr = read(hdf5,prop)[:,1:min(maxhits,nhits),:]
    return corr
end
function correlators(hdf5file,hdf5group,typeU,typeD,key;kws...)
    cU = correlators(hdf5file,hdf5group,typeU,key;kws...)
    cD = correlators(hdf5file,hdf5group,typeD,key;kws...)
    return @. (cU + cD)/2
end
