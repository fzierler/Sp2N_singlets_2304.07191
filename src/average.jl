function average_correlator(file::AbstractString,key::AbstractString,type::AbstractString;therm=0,step=1,kws...)
    corrs = correlators(file,type,key)[:,therm+1:step:end]
    average_correlator(key,corrs;kws...)
end
function average_correlator(hdf5file,hdf5group,key::AbstractString,type::AbstractString;therm=0,step=1,kws...)
    corrs = correlators(hdf5file,hdf5group,type,key)[:,therm+1:step:end]
    average_correlator(key,corrs;kws...)
end
average_correlator(key::String,corrs;autocor=true) = average_correlator(corrs;autocor=autocor)
function average_correlator(cor;autocor=true)
    N  = size(cor)[2]
    Nt = size(cor)[1]
    corr_mean = reshape(mean(cor,dims=2),Nt)
    # variance of the mean
    corr_var = reshape(var(cor,dims=2),Nt) ./ N
    # esitmate of covariance matrix of the mean
    cov_mat = cov(cor,dims=2) ./ N
    if autocor
        for i in 1:Nt
            τexp = autocorrelation_time(cor[i,:])
            corr_var[i] *= τexp
            cov_mat[i,i] *= τexp
        end
    end
    return corr_mean, corr_var, cov_mat
end
function average_plaquette(file;therm=0,autocor=true,step=1)
    plaqs = plaquettes(file;therm=therm,step=step)
    return average_plaquette(plaqs;autocor)
end
function average_plaquette(hdf5file,hdf5group;therm=0,autocor=true,step=1)
    plaqs = plaquettes(hdf5file,hdf5group;therm=therm,step=step)
    return average_plaquette(plaqs;autocor)
end
function average_plaquette(plaqs::Vector{S};autocor=true) where S <: Number
    p = mean(plaqs)
    Δp = std(plaqs)/sqrt(length(plaqs))
    if autocor
        τexp = autocorrelation_time(plaqs)
        Δp *= sqrt(τexp)
    end
    return p, Δp
end
function averagevectors!(corrs,T)
    tmp = zeros(T)
    keys_corr = keys(corrs[1])
    for dict in corrs
        for key in ["","g0","g5","g0g5"], suffix in ["","_re","_im"]
            k1 = key*"g1"*suffix
            k2 = key*"g2"*suffix
            k3 = key*"g3"*suffix
            if k1 in keys_corr && k2 in keys_corr && k3 in keys_corr
                @. tmp = (dict[k1] + dict[k2] + dict[k3])/3
                dict[k1] .= dict[k2] .= dict[k3] .= tmp
            end
        end
    end
end
#function average_correlator(file::String,key,typeU,typeD;kws...)
#    corrsU = correlators(file,typeU,key)
#    corrsD = correlators(file,typeD,key)
#    @. corrsU = (corrsU + corrsD)/2
#    average_correlator(key,corrsU)
#end
function average_correlator(hdf5file,hdf5group,key,typeU,typeD;kws...)
    corrsU = correlators(hdf5file,hdf5group,typeU,key)
    corrsD = correlators(hdf5file,hdf5group,typeD,key)
    @. corrsU = (corrsU + corrsD)/2
    average_correlator(key,corrsU)
end
