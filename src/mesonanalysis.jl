fitdecay(fit) = abs(getindex(fit.param,2))
function decayconstantmasses(fit1,fit2,fit3,error;constant=false)
    m = fitmass(fit1)
    f = fitdecay(fit1)
    c = constant ? fit1.param[3] : Float64(0)
    if error == :std
        # approximation by linearization around best fit
        Δm = stderror(fit1)[1]
        Δf = stderror(fit1)[2]
        Δc = constant ? stderror(fit1)[3] : Float64(0)
        return m, Δm, f, Δf, c, Δc
    elseif error == :bars
        # fitting the error bars
        Δm = abs(fitmass(fit3)-fitmass(fit2))/2
        Δf = abs(fitdecay(fit3)-fitdecay(fit2))/2
        Δc = constant ? abs(fit3.param[3]-fit2.param[3])/2 : Float64(0)
        return m, Δm, f, Δf, c, Δc
    else
        # pick larger error
        Δm1 = stderror(fit1)[1]
        Δf1 = stderror(fit1)[2]
        Δm2 = abs(fitmass(fit3)-fitmass(fit2))/2
        Δf2 = abs(fitdecay(fit3)-fitdecay(fit2))/2
        Δm = max(Δm1,Δm2)
        Δf = max(Δf1,Δf2)
        if constant
            Δc1 = stderror(fit1)[3]
            Δc2 = abs(fit3.param[3]-fit2.param[3])/2
            Δc = max(Δc1,Δc2)
        else
            Δc = Float64(0)
        end
        return mπ, Δmπ, fπ, Δfπ, c, Δc
    end
end
function decayconstantmasses(fit1;constant)
    m  = fitmass(fit1)
    f  = fitdecay(fit1)
    Δm = stderror(fit1)[1]
    Δf = stderror(fit1)[2]
    c  = constant ? fit1.param[3]     : Float64(0)
    Δc = constant ? stderror(fit1)[3] : Float64(0)
    return m, Δm, f, Δf, c, Δc
end
function _rescale_corrs!(corr,cov_v,cov_m,L)
    n1 = L^6/2 # from the norm used in HiReo
    n2 = L^3   # only keep a norm for the spatial volume
    @. corr  *= (n1/n2)
    @. cov_v *= (n1/n2)^2
    @. cov_m *= (n1/n2)^2
end
function _rescale_corrs!(corr,cov_v,L)
    n1 = L^6/2 # from the norm used in HiReo
    n2 = L^3   # only keep a norm for the spatial volume
    @. corr  *= (n1/n2)
    @. cov_v *= (n1/n2)^2
end
function _rescale_corrs!(corr,L)
    n1 = L^6/2 # from the norm used in HiReo
    n2 = L^3   # only keep a norm for the spatial volume
    @. corr  *= (n1/n2)
end
# DETERMINE MASSES AND DECAY CONSTANTS
#function meson_mass_decay(file::AbstractString,key,types...;constant=false,sign=+1,ncut,error,nexp2,kws...)
#    corr,cov_v,cov_m = average_correlator(file,key,types...;kws...)
#    L = latticesize(file)[2]
#    _rescale_corrs!(corr,cov_v,cov_m,L)
#    meson_mass_decay(corr,cov_v;ncut,sign,error,nexp2,constant)
#end
function meson_mass_decay_jackknife(file::AbstractString,key,type;kws...)
    corrs = correlators(file,type,key)
    L = latticesize(file)[2]
    _rescale_corrs!(corrs,L)
    meson_mass_decay_jackknife(corrs;kws...)
end
function meson_mass_decay_bootstrap(file::AbstractString,key,type;kws...)
    corrs = correlators(file,type,key)
    L = latticesize(file)[2]
    _rescale_corrs!(corrs,L)
    meson_mass_decay_bootstrap(corrs;nsamples=10_000,kws...)
end
function meson_mass_decay(hdf5file,hdf5group,key,types...;constant=false,sign=+1,ncut,error,nexp2,kws...)
    corr,cov_v,cov_m = average_correlator(hdf5file,hdf5group,key,types...;kws...)
    L = latticesize(hdf5file,hdf5group)[2]
    _rescale_corrs!(corr,cov_v,cov_m,L)
    meson_mass_decay(corr,cov_v;ncut,sign,error,nexp2,constant)
end
function meson_mass_decay_jackknife(hdf5file,hdf5group,key,type;kws...)
    corrs = correlators(hdf5file,hdf5group,type,key)
    L = latticesize(hdf5file,hdf5group)[2]
    _rescale_corrs!(corrs,L)
    meson_mass_decay_jackknife(corrs;kws...)
end
function meson_mass_decay_bootstrap(hdf5file,hdf5group,key,type;kws...)
    corrs = correlators(hdf5file,hdf5group,type,key)
    L = latticesize(hdf5file,hdf5group)[2]
    _rescale_corrs!(corrs,L)
    meson_mass_decay_bootstrap(corrs;nsamples=10_000,kws...)
end
function meson_mass_decay(corr,cov_v;sign=+1,constant=false,ncut,error,nexp2)
    if error == :bars
        fit1, fit2, fit3, model = fit_corr_bars(corr,cov_v,ncut;nexp2,sign,constant)
        mπ, Δmπ, fπ, Δfπ, c, Δc = decayconstantmasses(fit1,fit2,fit3,error;constant)
    elseif error == :none
        fit1, model = fit_corr(corr,cov_v,ncut;nexp2,sign,constant)
        mπ = fitmass(fit1)
        fπ = fitdecay(fit1)
        c = constant ? fit1.param[3] : Float64(0)
        return mπ, Float64(0), fπ, Float64(0), c, Float64(0)
    else
        fit1, model = fit_corr(corr,cov_v,ncut;nexp2,sign,constant)
        mπ, Δmπ, fπ, Δfπ, c, Δc = decayconstantmasses(fit1;constant)
    end
    return mπ, Δmπ, fπ, Δfπ, c, Δc
end
function meson_mass_decay_jackknife(corrs;sign=+1,constant=false,ncut,autocor,nexp2,kws...)
    T, N = size(corrs)
    # create arrays for decay constant
    corrs_delete1 = zeros(T,N-1)
    masses = zeros(N)
    decayc = zeros(N)
    constants = zeros(N)
    # set up jack-knife (one deletion)
    for i in 1:N
        for t in 1:T
            for j in 1:N
               (j < i) && (corrs_delete1[t,j]   = corrs[t,j])
               (j > i) && (corrs_delete1[t,j-1] = corrs[t,j])
            end
        end
        # perform averaging for fitting weights
        c_mean = reshape(mean(corrs_delete1,dims=2),T)
        cov_var = reshape(var(corrs_delete1,dims=2),T) ./ (N-1)
        #c_mean, cov_var = average_correlator(corrs_delete1;autocor)[1:2]
        fit = fit_corr(c_mean, cov_var,ncut;nexp2,sign,constant)[1]
        masses[i] = fitmass(fit)
        decayc[i] = fitdecay(fit)
        if constant
            constants[i] = fit.param[3]
        end
    end
    # now compute jack-knife best estimate and variance
    m, Δm  = mean(masses), sqrt(N-1)*std(masses,corrected=false)
    f, Δf  = mean(decayc), sqrt(N-1)*std(decayc,corrected=false)
    c, Δc  = mean(constants), sqrt(N-1)*std(constants,corrected=false)
    return m,Δm,f,Δf,c,Δc
end
function meson_mass_decay_bootstrap(corrs;nsamples=size(corrs)[2],sign=+1,constant=false,ncut,nexp2,kws...)
    T, N = size(corrs)
    # create arrays for decay constant
    corrs_sample = zeros(T,N)
    masses = zeros(nsamples)
    decayc = zeros(nsamples)
    constants = zeros(N)
    # set up jack-knife (one deletion)
    for i in 1:nsamples
        for j in 1:N
            ind = rand(1:N)
            corrs_sample[:,j] = corrs[:,ind]
        end
        # perform averaging for fitting weights
        c_mean = reshape(mean(corrs_sample,dims=2),T)
        cov_var = reshape(var(corrs_sample,dims=2),T)
        fit = fit_corr(c_mean, cov_var,ncut;nexp2,sign,constant)[1]
        masses[i] = fitmass(fit)
        decayc[i] = fitdecay(fit)
        if constant
            constants[i] = fit.param[3]
        end
    end
    # now compute jack-knife best estimate and variance
    m, Δm  = mean(masses), std(masses,corrected=false)
    f, Δf  = mean(decayc), std(decayc,corrected=false)
    c, Δc  = mean(constants), std(constants,corrected=false)
    return m,Δm,f,Δf,c,Δc
end
function meson_mass_decay_select(args...;error,kws...)
    if error == :hist || error == :histogram
        decay_mass_histogram(args...;kws...)
    elseif error == :jack || error == :jackknife
        meson_mass_decay_jackknife(args...;kws...)
    elseif error == :boot || error == :bootstrap
        println("using bootstrap")
        meson_mass_decay_bootstrap(args...;kws...)
    else
        meson_mass_decay(args...;error,kws...)
    end
end
#function groundstate_correlator(file::AbstractString,args...;sign=+1,kws...)
#    T = latticesize(file)[1]
#    t = collect(0:T-1)
#    m, Δm, f, Δf = meson_mass_decay_select(file,args...;sign,kws...)
#    model(t,m,f,T) = model_1st(t,T,m,f;sign)
#    model_err(t,m,f,T,Δf,Δm) = sqrt((model(t,m,1,T)*f*Δf)^2 + ((f^2/2*(exp(-m*t)+exp(-m*(T-t))) + f^2*m/2*(exp(-m*t)*t+exp(-m*(T-t))*(T-t)))*Δm)^2)
#    CI, ΔCI = model.(t,m,f,T), model_err.(t,m,f,T,Δf,Δm)
#    return CI, ΔCI
#end
function groundstate_correlator(hdf5file::AbstractString,hdf5group::AbstractString,args...;sign=+1,kws...)
    T = latticesize(hdf5file,hdf5group)[1]
    t = collect(0:T-1)
    m, Δm, f, Δf = meson_mass_decay_select(hdf5file,hdf5group,args...;sign,kws...)
    model(t,m,f,T) = model_1st(t,T,m,f;sign)
    model_err(t,m,f,T,Δf,Δm) = sqrt((model(t,m,1,T)*f*Δf)^2 + ((f^2/2*(exp(-m*t)+exp(-m*(T-t))) + f^2*m/2*(exp(-m*t)*t+exp(-m*(T-t))*(T-t)))*Δm)^2)
    CI, ΔCI = model.(t,m,f,T), model_err.(t,m,f,T,Δf,Δm)
    return CI, ΔCI
end
function groundstate_correlator(corr,cut;sign=+1,kws...)
    CC, CCVar = average_correlator(corr;kws...)[1:2]
    m, Δm, f, Δf = meson_mass_decay(CC,CCVar;ncut=cut,sign,error=:std,nexp2=false)
    T = length(CC)
    t = collect(0:T-1)
    model(t,m,f,T) = model_1st(t,T,m,f;sign)
    model_err(t,m,f,T,Δf,Δm) = sqrt((model(t,m,1,T)*f*Δf)^2 + ((f^2/2*(exp(-m*t)+exp(-m*(T-t))) + f^2*m/2*(exp(-m*t)*t+exp(-m*(T-t))*(T-t)))*Δm)^2)
    CI, ΔCI = model.(t,m,f,T), model_err.(t,m,f,T,Δf,Δm)
    return CI, ΔCI
end
# jackknife related methods for singlets
function correlator_deriv(C,ΔC)
    T = length(C)
    Cd  = zero(C)
    ΔCd = zero(C)
    for i in 1:T
        Cd[i]  = ( C[mod1(i-1,T)] -  C[mod1(i+1,T)])/2
        ΔCd[i] = sqrt((ΔC[mod1(i-1,T)]^2 + ΔC[mod1(i+1,T)]^2)/4)
    end
    return Cd, ΔCd
end
function apply_jackknife(obs;ignoreNaN=false)
    N  = size(obs)[1]
    if ignoreNaN
        O  = vec(nanmean(obs,dims=1))
        ΔO = vec(sqrt(N-1)*nanstd(obs,corrected=false,dims=1))
    else
        O  = vec(mean(obs,dims=1))
        ΔO = vec(sqrt(N-1)*std(obs,corrected=false,dims=1))
    end
    return O, ΔO
end
function apply_bootstrap(obs;ignoreNaN=false)
    if ignoreNaN
        O  = vec(nanmean(obs,dims=1))
        ΔO = vec(nanstd(obs,corrected=false,dims=1))
    else
        O  = vec(mean(obs,dims=1))
        ΔO = vec(std(obs,corrected=false,dims=1))
    end
    return O, ΔO
end
function singlet_jackknife(C_con_MC,C_dis_MC,cut,fitint;deriv=true,gs_sub=true,sigma=false,constant=false)
    # set-up jackknife
    N1, T = size(C_dis_MC)
    T, N2 = size(C_con_MC)
    N = min(N1,N2)
    mη_jk   = zeros(N)
    Cd_jk   = zeros(N,T÷1)
    C_jk    = zeros(N,T÷1)
    meff_jk = zeros(N,T÷2)
    const_η_jk  = zeros(N)
    C_con_MC_JK = zeros(T,N-1)
    C_dis_MC_JK = zeros(N-1,T)
    for i in 1:N
        for t in 1:T
            for j in 1:N
            (j < i) && (C_con_MC_JK[t,j]   = C_con_MC[t,j])
            (j < i) && (C_dis_MC_JK[j,t]   = C_dis_MC[j,t])
            (j > i) && (C_con_MC_JK[t,j-1] = C_con_MC[t,j])
            (j > i) && (C_dis_MC_JK[j-1,t] = C_dis_MC[j,t])
            end
        end
        # Average connected and disconnect pieces
        τ = 1 #autocorrelation_time(plaquettes(fileS);correct=true)
        C  = mean(C_dis_MC_JK,dims=1)[1,:]
        ΔC = sqrt.(var(C_dis_MC_JK,dims=1)[1,:]/((N-1)/τ))
        # full connected diagrams
        Cc = mean(C_con_MC_JK,dims=2)[:,1]
        ΔCc= sqrt.(var(C_con_MC_JK,dims=2)[:,1]/((N-1)/τ))
        # obtain ground state correlator of connected part
        if gs_sub
            CI, ΔCI = groundstate_correlator(C_con_MC_JK,cut;sign=+1,autocor=false)
        end
        # combine the diagramsgs_sub
        Nf=2
        sign = sigma ? +1 : -1
        if gs_sub
            improved  = @. CI + sign*Nf*C
            Δimproved = @. sqrt(ΔCI^2 + Nf^2*ΔC^2)
        else
            improved  = @. Cc + sign*Nf*C
            Δimproved = @. sqrt(ΔCc^2 + Nf^2*ΔC^2)
        end
        # derivatiuve method
        Cd, ΔCd = correlator_deriv(improved,Δimproved)
        if deriv
            mη, Δmη, fη, Δfη, cη, Δcη = meson_mass_decay(Cd,ΔCd.^2;ncut=fitint,sign=-1,error=:none,nexp2=false,constant)
            meffη = implicit_meff(Cd,sign=-1)
        else
            try
                if constant
                    fit_constant = fitint #(T÷2-4,T÷2)
                    mη, Δmη, fη, Δfη, cη, Δcη = meson_mass_decay(improved,Δimproved.^2;ncut=fit_constant,sign=+1,error=:std,nexp2=false,constant=true)
                    improved = improved .- cη
                    mη, Δmη, fη, Δfη = meson_mass_decay(improved,Δimproved.^2;ncut=fitint,sign=+1,error=:std,nexp2=false,constant=false)[1:4]
                else
                    mη, Δmη, fη, Δfη, cη, Δcη = meson_mass_decay(improved,Δimproved.^2;ncut=fitint,sign=+1,error=:std,nexp2=false,constant)
                end
            catch
                mη, Δmη, fη, Δfη, cη, Δcη = 0, 0, 0, 0, 0 ,0
            end
            if constant
                meffη = implicit_meff(improved .- cη ,sign=+1)
            else
                meffη = implicit_meff(improved,sign=+1)
            end
        end
        # fill jk samples
        Cd_jk[i,:]  .= Cd[:]
        C_jk[i,:]   .= improved[:]
        meff_jk[i,:].= meffη[:]
        mη_jk[i]     = mη
        const_η_jk[i]= cη
    end
    mη, Δmη = first.(apply_jackknife(mη_jk))
    # Ignoring NaNs is only used for the effective masses
    meffη, Δmeffη = apply_jackknife(meff_jk;ignoreNaN=true)
    Cd, ΔCd = apply_jackknife(Cd_jk)
    C, ΔC = apply_jackknife(C_jk)
    const_η, Δconst_η = first.(apply_jackknife(const_η_jk))
    return mη, Δmη, Cd, ΔCd, C, ΔC, meffη, Δmeffη, const_η, Δconst_η
end
function singlet_jackknife_corr(C_con_MC,C_dis_MC,cut;gs_sub=true,sigma=false)
    # set-up jackknife
    N1, T = size(C_dis_MC)
    T, N2 = size(C_con_MC)
    N = min(N1,N2)
    Cd_jk   = zeros(N,T÷1)
    C_jk    = zeros(N,T÷1)
    C_con_MC_JK = zeros(T,N-1)
    C_dis_MC_JK = zeros(N-1,T)
    for i in 1:N
        for t in 1:T
            for j in 1:N
            (j < i) && (C_con_MC_JK[t,j]   = C_con_MC[t,j])
            (j < i) && (C_dis_MC_JK[j,t]   = C_dis_MC[j,t])
            (j > i) && (C_con_MC_JK[t,j-1] = C_con_MC[t,j])
            (j > i) && (C_dis_MC_JK[j-1,t] = C_dis_MC[j,t])
            end
        end
        # Average connected and disconnect pieces
        τ = 1 #autocorrelation_time(plaquettes(fileS);correct=true)
        C  = mean(C_dis_MC_JK,dims=1)[1,:]
        ΔC = sqrt.(var(C_dis_MC_JK,dims=1)[1,:]/((N-1)/τ))
        # full connected diagrams
        Cc = mean(C_con_MC_JK,dims=2)[:,1]
        ΔCc= sqrt.(var(C_con_MC_JK,dims=2)[:,1]/((N-1)/τ))
        # obtain ground state correlator of connected part
        if gs_sub
            CI, ΔCI = groundstate_correlator(C_con_MC_JK,cut;sign=+1,autocor=false)
        end
        # combine the diagramsgs_sub
        Nf=2
        sign = sigma ? +1 : -1
        if gs_sub
            improved  = @. CI + sign*Nf*C
            Δimproved = @. sqrt(ΔCI^2 + Nf^2*ΔC^2)
        else
            improved  = @. Cc + sign*Nf*C
            Δimproved = @. sqrt(ΔCc^2 + Nf^2*ΔC^2)
        end
        # derivatiuve method
        Cd, ΔCd = correlator_deriv(improved,Δimproved)
        # fill jk samples
        Cd_jk[i,:]  .= Cd[:]
        C_jk[i,:]   .= improved[:]
    end
    Cd, ΔCd = apply_jackknife(Cd_jk)
    C, ΔC = apply_jackknife(C_jk)
    return C, ΔC, Cd, ΔCd
end

function singlet_bootstrap(C_con_MC,C_dis_MC,cut,fitint;deriv=true,gs_sub=true,sigma=false,constant=false)
    # set-up jackknife
    N1, T = size(C_dis_MC)
    T, N2 = size(C_con_MC)
    N = min(N1,N2)
    nsamples = 10N
    mη_BS   = zeros(nsamples)
    Cd_BS   = zeros(nsamples,T÷1)
    C_BS    = zeros(nsamples,T÷1)
    meff_BS = zeros(nsamples,T÷2)
    const_η_BS  = zeros(nsamples)
    C_con_MC_BS = zeros(T,N)
    C_dis_MC_BS = zeros(N,T)
    for i in 1:nsamples
        for j in 1:N
            ind = rand(1:N)
            C_con_MC_BS[:,j] = C_con_MC[:,ind]
            C_dis_MC_BS[j,:] = C_dis_MC[ind,:]
        end
        # Average connected and disconnect pieces
        τ = 1 #autocorrelation_time(plaquettes(fileS);correct=true)
        C  = mean(C_dis_MC_BS,dims=1)[1,:]
        ΔC = sqrt.(var(C_dis_MC_BS,dims=1)[1,:]/(N/τ))
        # full connected diagrams
        Cc = mean(C_con_MC_BS,dims=2)[:,1]
        ΔCc= sqrt.(var(C_con_MC_BS,dims=2)[:,1]/(N/τ))
        # obtain ground state correlator of connected part
        if gs_sub
            CI, ΔCI = groundstate_correlator(C_con_MC_BS,cut;sign=+1,autocor=false)
        end
        # combine the diagramsgs_sub
        Nf=2
        sign = sigma ? +1 : -1
        if gs_sub
            improved  = @. CI + sign*Nf*C
            Δimproved = @. sqrt(ΔCI^2 + Nf^2*ΔC^2)
        else
            improved  = @. Cc + sign*Nf*C
            Δimproved = @. sqrt(ΔCc^2 + Nf^2*ΔC^2)
        end
        # derivatiuve method
        Cd, ΔCd = correlator_deriv(improved,Δimproved)
        if deriv
            mη, Δmη, fη, Δfη, cη, Δcη = meson_mass_decay(Cd,ΔCd.^2;ncut=fitint,sign=-1,error=:std,nexp2=false,constant)
            meffη = implicit_meff(Cd,sign=-1)
        else
            mη, Δmη, fη, Δfη, cη, Δcη = meson_mass_decay(improved,Δimproved.^2;ncut=fitint,sign=+1,error=:std,nexp2=false,constant)
            if constant
                meffη = implicit_meff(improved .- cη ,sign=+1)
            else
                meffη = implicit_meff(improved,sign=+1)
            end
        end
        # fill jk samples
        Cd_BS[i,:]  .= Cd[:]
        C_BS[i,:]   .= improved[:]
        meff_BS[i,:].= meffη[:]
        mη_BS[i]     = mη
        const_η_BS[i]= cη
    end
    return Cd_BS, C_BS, meff_BS, mη_BS, const_η_BS
end
