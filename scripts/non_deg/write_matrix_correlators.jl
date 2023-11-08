using HDF5
using DelimitedFiles
using HiRepAnalysis

function fermion_masses_from_filename(file)
    p1 = last(findfirst("m1",file))
    p2 = first(findfirst("m2",file))
    p3 = last(findfirst("m2",file))
    # create array of masses for matching of output
    m  = [file[p1+1:p2-1],file[p3+1:end]]
    return m
end
function unbiased_estimator(discon1,discon2;nhits,rescale=1)
    # (1) average over different hits
    # (2) average over all time separations
    # (3) normalize wrt. time and hit average
    @assert size(discon1)[2] == size(discon2)[2] == nhits
    nconf, nhits, T = size(discon1)
    timavg = zeros(eltype(discon1),(nconf,T))
    norm   = T*nhits^2
    for t in 1:T
        for t0 in 1:T
            Δt = mod(t-t0,T)
            @inbounds for hit1 in 1:nhits, hit2 in 1:nhits
                for conf in 1:nconf
                    timavg[conf,Δt+1] += discon1[conf,hit1,t]*discon2[conf,hit2,t0]
                end
            end
        end
    end
    @. timavg = rescale*timavg/norm
    # transpose the matrix so that it has the same layout as the connected pieces
    return permutedims(timavg)
end
function unbiased_estimator(discon;nhits,rescale=1)
    # (1) average over different hits
    # (2) average over all time separations
    # (3) normalize wrt. time and hit average
    @assert size(discon)[2] == nhits
    nconf, nhits, T = size(discon)
    timavg = zeros(eltype(discon),(nconf,T))
    norm   = T*div(nhits,2)^2
    hitsd2 = div(nhits,2)
    for t in 1:T
        for t0 in 1:T
            Δt = mod(t-t0,T)
            @inbounds for hit1 in 1:hitsd2, hit2 in hitsd2+1:nhits
                for conf in 1:nconf
                    timavg[conf,Δt+1] += discon[conf,hit1,t]*discon[conf,hit2,t0]
                end
            end
        end
    end
    @. timavg = rescale*timavg/norm
    return permutedims(timavg)
end
function rescale_connected!(corr,L)
    n1 = L^6/2 # from the norm used in HiReo
    n2 = L^3   # only keep a norm for the spatial volume
    @. corr  *= (n1/n2)
end
function read_hdf5_zenodo(ens,Γ)
    # parse fermion masses from filename
    logfile = splitpath(read(ens,"out_spectrum_discon/logfile"))[4]
    m1, m2 = fermion_masses_from_filename(logfile)

    connU = read(ens,"out_spectrum_with_pcac/SEMWALL_D TRIPLET_$Γ")
    connD = read(ens,"out_spectrum_with_pcac/SEMWALL_U TRIPLET_$Γ")
    discU = read(ens,"out_spectrum_discon/DISCON_SEMWALL SINGLET_$(Γ)_disc_re_$m1")
    discD = read(ens,"out_spectrum_discon/DISCON_SEMWALL SINGLET_$(Γ)_disc_re_$m2")
    # as a cross check get the number of stochastic sources from the hdf5 file
    nhits = read(ens,"out_spectrum_discon/sources")
    nconf = length(read(ens,"out_spectrum_discon/configurations"))
    T,L = read(ens,"out_spectrum_discon/lattice")[1:2]

    #rescale disconnected pieces to match the common normalisation
    rescale_disc = (L^3)^2 /L^3
    # use the stochastic estimator of arXiv:1607.06654 eq. (14)
    discUU = unbiased_estimator(discU;nhits,rescale=rescale_disc)
    discDD = unbiased_estimator(discD;nhits,rescale=rescale_disc)
    discUD = unbiased_estimator(discU,discD;nhits,rescale=rescale_disc)
    # rescale now the connected pieces appropriately
    rescale_connected!(connU,L)
    rescale_connected!(connD,L)
    return connU, connD, discUU, discDD, discUD
end
# choose two sets of operators which shoud lead to identical results. 
# CASE A: η': (uΓu + dΓd)/√2 
#         π⁰: (uΓu - dΓd)/√2
# CASE B: O1: uΓu 
#         O2: dΓd  

# Operators: η': (uΓu + dΓd)/√2 
#            π⁰: (uΓu - dΓd)/√2
function correlation_matrix_basisA(connU,connD,discUU,discDD,discUD)
    corr_matrix = zeros(eltype(connU),(2,2,size(connU)...))
    # create the correlators and cross-correlators of CASE 1
    @. corr_matrix[1,1,:,:] = (connU + connD)/2 - discUD - (discUU + discDD)/2
    @. corr_matrix[2,2,:,:] = (connU + connD)/2 + discUD - (discUU + discDD)/2
    @. corr_matrix[1,2,:,:] = -(connD - connU)/2 - (discUU - discDD)/2
    @. corr_matrix[2,1,:,:] = -(connD - connU)/2 - (discUU - discDD)/2
    return corr_matrix 
end
# Operators: O1: uΓu 
#            O2: dΓd  
function correlation_matrix_basisB(connU,connD,discUU,discDD,discUD)
    corr_matrix = zeros(eltype(connU),(2,2,size(connU)...))
    # create the correlators and cross-correlators of CASE 2
    @. corr_matrix[1,1,:,:] = -connU + discUU
    @. corr_matrix[2,2,:,:] = -connD + discDD
    @. corr_matrix[1,2,:,:] = discUD
    @. corr_matrix[2,1,:,:] = discUD
    return corr_matrix 
end
function correlator_derivative(c::AbstractArray;t_dim=1)
    c_deriv = similar(c)
    # Get the number of Euclidean timeslices in the specified dimension
    T = size(c)[t_dim]
    # get the overall number of dimensions for contructing a unit vector in the
    # temporal direction.
    n = ndims(c)
    δt(i) = ifelse(i == t_dim,1,0)
    # ntuple applies the function in the first argument δt to every index `1:n`
    # Thus, this corresponds to a delta-step in the temporal direction, i.e. a 
    # unit vector in Euclidean time. 
    unit_t = CartesianIndex(ntuple(δt,n))
    # Use julia's Cartesian indices to index arrays of arbitrary dimension 
    for i in CartesianIndices(c)
        # special case the first and the last index: 
        if i[t_dim] == 1
           c_deriv[i] = c[i+unit_t] - c[i]
        elseif i[t_dim] == T
           c_deriv[i] = c[i] - c[i-unit_t]
        else
           c_deriv[i] = (c[i+unit_t] - c[i-unit_t])/2
        end    
    end
    return c_deriv
end
function write_matrix_correlators(hdf5dir,filename,ens)
    ispath(hdf5dir) || mkpath(hdf5dir)
    file = h5open(joinpath(hdf5dir,filename),"w")

    all_keys = keys(ens["out_spectrum_discon"])
    keys_without_correlators = filter(!contains("SINGLET"),all_keys)

    # save data that has not been changed
    for k in keys_without_correlators
        contains(k, "quarkmasses") && break
        file[k] = read(ens["out_spectrum_discon"],k)
    end
    # write masses to file
    m = fermion_masses_from_filename(splitext(filename)[1])
    file["quarkmasses"] = m

    # get connected and disconnected piece of pseuooscalar correlator
    connU, connD, discUU, discDD, discUD = read_hdf5_zenodo(ens,"g5")

    # set up correlator matrices
    corr1 = correlation_matrix_basisA(connU,connD,discUU,discDD,discUD)
    corr2 = correlation_matrix_basisB(connU,connD,discUU,discDD,discUD)
    corr1_deriv = correlator_derivative(corr1;t_dim=3) 
    corr2_deriv = correlator_derivative(corr1;t_dim=3) 

    write_matrix_julia(file,corr1,corr2,corr1_deriv,corr2_deriv)
    close(file)
end
function write_matrix_julia(file,corr1,corr2,corr1_deriv,corr2_deriv)
    file["correlator_matrix_A"] = corr1
    file["correlator_matrix_B"] = corr2
    file["correlator_matrix_A_deriv"] = corr1_deriv
    file["correlator_matrix_B_deriv"] = corr2_deriv
end
function write_all_matrix_correlators(dat_file,prm_file;hdf5dir,kws...)
    # open HDF5 file that contains all data from the first singlet paper
    file_id = h5open(dat_file)
    parameters = readdlm(prm_file,';',skipstart=1)

    # restrict ourselves to the hdf5 ensembles with non-degenerate fermion masses.
    run_labels = parameters[:,1]
    Sp4runs = getindex.(Ref(file_id),run_labels)

    for (ind,h5group) in enumerate(Sp4runs)
        filename = basename(parameters[ind,1])*".h5"
        m = fermion_masses_from_filename(basename(parameters[ind,1]))
        if m[1] != m[2] 
            write_matrix_correlators(hdf5dir,filename,h5group;kws...)
        end
    end
end

dat_file = "input_data/data.hdf5"
prm_file = "input/parameters/param_non_deg.csv"
out_file = "output/matrix_correlators"
write_all_matrix_correlators(dat_file,prm_file;hdf5dir=out_file)