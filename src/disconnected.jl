# This is for ouput files from /Spectrum/
function flatten_disc(d::Vector{Vector{S}},nhits;rescale=1) where S
    # the rounding automatically removes incomplete calcluations
    n = div(length(d),nhits)
    T = length(d[1])
    m = zeros(eltype(d[1]),(n, nhits, T))
    for i in 1:n, j in 1:nhits, t in 1:T
        m[i,j,t] = d[(i-1)*nhits + j][t]
    end
    @. m = rescale*m
    return m
end
function flatten_disc(d::Matrix{S},nhits;rescale=1) where S
    # the rounding automatically removes incomplete calcluations
    T, N = size(d)
    n = div(N,nhits)
    m = zeros(eltype(d[1]),(n, nhits, T))
    for i in 1:n, j in 1:nhits, t in 1:T
        m[i,j,t] = d[t,(i-1)*nhits + j]
    end
    @. m = rescale*m
    return m
end
# time averaged construction of disconnected diagrams
function hit_time_average_vacuum(m;rescale=1)
    nconf, nhits, T = size(m)
    timavg = zeros(eltype(m),(nconf,T))
    @inbounds for t in 1:T, hit1 in 1:nhits, conf in 1:nconf
        timavg[conf,t] += m[conf,hit1,t]
    end
    @. timavg = timavg*(rescale/nhits)
    return timavg
end
function hit_time_average_disconnected(m;rescale=1)
    # (1) average over different hits
    # (2) average over all time separations
    # (3) normalize wrt. time and hit average
    nconf, nhits, T = size(m)
    timavg = zeros(eltype(m),(nconf,T))
    norm   = T*div(nhits,2)^2
    hitsd2 = div(nhits,2)
    for t in 1:T
        for t0 in 1:T
            Δt = mod(t-t0,T)
            @inbounds for hit1 in 1:hitsd2, hit2 in hitsd2+1:nhits
                for conf in 1:nconf
                    timavg[conf,Δt+1] += conj(m[conf,hit1,t])*m[conf,hit2,t0]
                end
            end
        end
    end
    @. timavg = rescale*timavg/norm
    return timavg
end
function hit_time_average_disconnected_vacuum(m;rescale=1)
    nconf, nhits, T = size(m)
    hitsd2 = div(nhits,2)
    # (1) average over different hits
    # (2) average over all time separations
    # (3) normalize wrt. time and hit average
    v1 = hit_time_average_vacuum(m[:,1:hitsd2,:])
    v1 = dropdims(mean(v1,dims=1),dims=1)
    v1 = mean(v1[2:end])
    #@show v1^2*rescale
    v2 = hit_time_average_vacuum(m[:,hitsd2+1:end,:])
    v2 = dropdims(mean(v2,dims=1),dims=1)
    v2 = mean(v2[2:end])
    # now fully disconnected
    timavg = zeros(eltype(m),(nconf,T))
    norm   = T*div(nhits,2)^2
    for t in 1:T
        for t0 in 1:T
            Δt = mod(t-t0,T)
            @inbounds for hit1 in 1:hitsd2, hit2 in hitsd2+1:nhits
                for conf in 1:nconf
                    #timavg[conf,Δt+1] += conj(m[conf,hit1,t]-v1[t])*(m[conf,hit2,t0]-v2[t0])
                    #timavg[conf,Δt+1] += conj(m[conf,hit1,t])*(m[conf,hit2,t0]) - v1*v2
                    timavg[conf,Δt+1] += conj(m[conf,hit1,t]-v1)*(m[conf,hit2,t0]-v2)
                end
            end
        end
    end
    @. timavg = rescale*timavg/norm
    return timavg
end
function hit_time_average_disconnected_vacuum(m1,m2;rescale=1)
    v1 = hit_time_average_vacuum(m1)
    v1 = dropdims(mean(v1,dims=1),dims=1)
    #v1 = mean(v1[2:end])
    v2 = hit_time_average_vacuum(m2)
    v2 = dropdims(mean(v2,dims=1),dims=1)
    # (1) average over different hits
    # (2) average over all time separations
    # (3) normalize wrt. time and hit average
    nconf, nhits, T = size(m1)
    timavg = zeros(eltype(m1),(nconf,T))
    norm   = T*nhits^2
    for t in 1:T
        for t0 in 1:T
            Δt = mod(t-t0,T)
            @inbounds for hit1 in 1:nhits, hit2 in 1:nhits
                for conf in 1:nconf
                    timavg[conf,Δt+1] += conj(m1[conf,hit1,t]-v1[t])*(m2[conf,hit2,t0]-v2[t])
                end
            end
        end
    end
    @. timavg = rescale*timavg/norm
    return timavg
end
function hit_time_average_disconnected(m1,m2;rescale=1)
    # (1) average over different hits
    # (2) average over all time separations
    # (3) normalize wrt. time and hit average
    nconf, nhits, T = size(m1)
    timavg = zeros(eltype(m1),(nconf,T))
    norm   = T*nhits^2
    for t in 1:T
        for t0 in 1:T
            Δt = mod(t-t0,T)
            @inbounds for hit1 in 1:nhits, hit2 in 1:nhits
                for conf in 1:nconf
                    timavg[conf,Δt+1] += conj(m1[conf,hit1,t])*m2[conf,hit2,t0]
                end
            end
        end
    end
    @. timavg = rescale*timavg/norm
    return timavg
end
# This is for ouput files from /Spectrum/
function disconnected_eta(file,type,hits;key="g5_disc_re",maxhits=hits,kws...)
    c = correlators(file,type,key,hits;maxhits,withsource=true,average=false,filterkey=false,key_pattern=key)
    τ = max(1.0,autocorrelation_time(plaquettes(file)))
    _disconnected_eta(c,τ;kws...)
end
function disconnected_eta(hdf5file,hdf5group,type,hits;key="g5_disc_re",maxhits=hits,kws...)
    c = correlators(hdf5file,hdf5group,type,key,hits;maxhits,withsource=true,average=false,filterkey=false,key_pattern=key)
    τ = max(1.0,autocorrelation_time(plaquettes(hdf5file,hdf5group)))
    _disconnected_eta(c,τ;kws...)
end
# Function w/o Monte-Carlo average for jackknife
function disconnected_eta_MC(file,type,hits;maxhits=hits,vsub=false,key="g5_disc_re",kws...)
    c = correlators(file,type,key,hits;withsource=true,average=false,filterkey=false,key_pattern=key,maxhits)
    return vsub ? hit_time_average_disconnected_vacuum(c;kws...) : hit_time_average_disconnected(c;kws...)
end
function disconnected_eta_MC(hdf5file,hdf5group,type,hits;maxhits=hits,vsub=false,key="g5_disc_re",kws...)
    c = correlators(hdf5file,hdf5group,type,key,hits;withsource=true,average=false,filterkey=false,key_pattern=key,maxhits)
    return vsub ? hit_time_average_disconnected_vacuum(c;kws...) : hit_time_average_disconnected(c;kws...)
end
function disconnected_eta(file,type,hits,m1,m2;key="g5_disc_re",maxhits=hits,therm=0,kws...)
    c1 = correlators(file,type,key,hits;withsource=true,average=false,masses=true,mass=m1,filterkey=false,key_pattern=key,maxhits)
    c2 = correlators(file,type,key,hits;withsource=true,average=false,masses=true,mass=m2,filterkey=false,key_pattern=key,maxhits)
    τ = max(1.0,autocorrelation_time(plaquettes(file)))
    _disconnected_nondeg(c1,c2,τ;sign=+1,kws...)
end
function disconnected_eta(hdf5file,hdf5group,type,hits,m1,m2;key="g5_disc_re",maxhits=hits,therm=0,kws...)
    c1 = correlators(hdf5file,hdf5group,type,key,hits;withsource=true,average=false,masses=true,mass=m1,filterkey=false,key_pattern=key,maxhits)
    c2 = correlators(hdf5file,hdf5group,type,key,hits;withsource=true,average=false,masses=true,mass=m2,filterkey=false,key_pattern=key,maxhits)
    τ = max(1.0,autocorrelation_time(plaquettes(hdf5file,hdf5group)))
    _disconnected_nondeg(c1,c2,τ;sign=+1,kws...)
end
function disconnected_pi0(file,type,hits,m1,m2;key="g5_disc_re",maxhits=hits,kws...)
    c1 = correlators(file,type,key,hits;withsource=true,average=false,masses=true,mass=m1,filterkey=false,key_pattern=key,maxhits)
    c2 = correlators(file,type,key,hits;withsource=true,average=false,masses=true,mass=m2,filterkey=false,key_pattern=key,maxhits)
    τ = max(1.0,autocorrelation_time(plaquettes(file)))
    _disconnected_nondeg(c1,c2,τ;sign=-1,kws...)
end
function disconnected_pi0(hdf5file,hdf5group,type,hits,m1,m2;key="g5_disc_re",maxhits=hits,kws...)
    c1 = correlators(hdf5file,hdf5group,type,key,hits;withsource=true,average=false,masses=true,mass=m1,filterkey=false,key_pattern=key,maxhits)
    c2 = correlators(hdf5file,hdf5group,type,key,hits;withsource=true,average=false,masses=true,mass=m2,filterkey=false,key_pattern=key,maxhits)
    τ = max(1.0,autocorrelation_time(plaquettes(hdf5file,hdf5group)))
    _disconnected_nondeg(c1,c2,τ;sign=-1,kws...)
end
function disconnected_eta_MC(file,type,hits,m1,m2;key="g5_disc_re",maxhits=hits,therm=0,kws...)
    c1 = correlators(file,type,key,hits;withsource=true,average=false,masses=true,mass=m1,filterkey=false,key_pattern=key,maxhits)
    c2 = correlators(file,type,key,hits;withsource=true,average=false,masses=true,mass=m2,filterkey=false,key_pattern=key,maxhits)
    τ = max(1.0,autocorrelation_time(plaquettes(file)))
    _disconnected_nondeg_MC(c1,c2,τ;sign=+1,kws...)
end
function disconnected_eta_MC(hdf5file,hdf5group,type,hits,m1,m2;key="g5_disc_re",maxhits=hits,therm=0,kws...)
    c1 = correlators(hdf5file,hdf5group,type,key,hits;withsource=true,average=false,masses=true,mass=m1,filterkey=false,key_pattern=key,maxhits)
    c2 = correlators(hdf5file,hdf5group,type,key,hits;withsource=true,average=false,masses=true,mass=m2,filterkey=false,key_pattern=key,maxhits)
    τ = max(1.0,autocorrelation_time(plaquettes(hdf5file,hdf5group)))
    _disconnected_nondeg_MC(c1,c2,τ;sign=+1,kws...)
end
function disconnected_pi0_MC(file,type,hits,m1,m2;key="g5_disc_re",maxhits=hits,kws...)
    c1 = correlators(file,type,key,hits;withsource=true,average=false,masses=true,mass=m1,filterkey=false,key_pattern=key,maxhits)
    c2 = correlators(file,type,key,hits;withsource=true,average=false,masses=true,mass=m2,filterkey=false,key_pattern=key,maxhits)
    τ = max(1.0,autocorrelation_time(plaquettes(file)))
    _disconnected_nondeg_MC(c1,c2,τ;sign=-1,kws...)
end
function disconnected_pi0_MC(hdf5file,hdf5group,type,hits,m1,m2;key="g5_disc_re",maxhits=hits,kws...)
    c1 = correlators(hdf5file,hdf5group,type,key,hits;withsource=true,average=false,masses=true,mass=m1,filterkey=false,key_pattern=key,maxhits)
    c2 = correlators(hdf5file,hdf5group,type,key,hits;withsource=true,average=false,masses=true,mass=m2,filterkey=false,key_pattern=key,maxhits)
    τ = max(1.0,autocorrelation_time(plaquettes(hdf5file,hdf5group)))
    _disconnected_nondeg_MC(c1,c2,τ;sign=-1,kws...)
end
# Methods for already parsed and pre-processed data
function _disconnected_eta(flattened,τ;maxconf=typemax(Int),vsub=false,kws...)
    timavg = vsub ? hit_time_average_disconnected_vacuum(flattened;kws...) : hit_time_average_disconnected(flattened;kws...)
    nconf  = size(timavg)[1]
    n = min(nconf,maxconf)
    timavg = timavg[1:n,:]
    C  = mean(timavg,dims=1)[1,:]
    ΔC = sqrt.(var(timavg,dims=1)[1,:]/n)
    return C, ΔC*sqrt(τ)
end
function _disconnected_nondeg(flat1,flat2,τ;sign,kws...)
    t = _disconnected_nondeg_MC(flat1,flat2,τ;sign,kws...)
    n = size(t)[1]
    C = mean(t[1:n,:],dims=1)[1,:]
    ΔC = sqrt.(var(t[1:n,:],dims=1)[1,:])/sqrt(n/τ)
    return C, ΔC
end
function _disconnected_nondeg_MC(flat1,flat2,τ;sign,vsub,kws...)
    if vsub
        t11 = hit_time_average_disconnected_vacuum(flat1;kws...)
        t22 = hit_time_average_disconnected_vacuum(flat2;kws...)
        t12 = hit_time_average_disconnected_vacuum(flat1,flat2;kws...)
    else
        t11 = hit_time_average_disconnected(flat1;kws...)
        t22 = hit_time_average_disconnected(flat2;kws...)
        t12 = hit_time_average_disconnected(flat1,flat2;kws...)
    end
    t = @. (t11 + t22 + sign*2t12)/4
    return t
end
# Only obtain vacuum constant
function disconnected_vacuum(file,type,hits;key="g5_disc_re",maxhits=hits,rescale)
    c = correlators(file,type,key,hits;maxhits,withsource=true,average=false,filterkey=false,key_pattern=key)
    v = hit_time_average_vacuum(c)
    v = mean(v)
    return v*rescale
end
function disconnected_vacuum(hdf5file,hdf5group,type,hits;key="g5_disc_re",maxhits=hits,rescale)
    c = correlators(hdf5file,hdf5group,type,key,hits;maxhits,withsource=true,average=false,filterkey=false,key_pattern=key)
    v = hit_time_average_vacuum(c)
    v = mean(v)
    return v*rescale
end
