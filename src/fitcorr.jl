model_1st(t,T,m,f;sign=+1,)  = (exp(-m*t)+sign*exp(-m*(T-t)))*f^2*m/2
model_2nd(t,T,m,f;sign=+1)  = (exp(-m*t)+sign*exp(-m*(T-t)))*f/2
fitmass(fit)=first(coef(fit))
function _aply_cut(ncut,T)
    ncut = isa(ncut,Int) ? (ncut,div(T,2)) : ncut
    @assert isa(ncut,Tuple)
    @assert length(ncut)==2
    return ncut[1]:ncut[2]
end
function fit_corr_1exp(c,w,ncut;sign=+1)
    T = length(c)
    t = collect(1:T)
    cut = _aply_cut(ncut,T)
    if ndims(w)==2
        weight = w[cut,cut]
        weight = (weight + weight')/2
    else
        weight = w[cut]
    end
    @. model(t,p) = model_1st(t-1,T,p[1],p[2];sign)
    fit = curve_fit(model,t[cut],c[cut],weight,ones(2))
    fit.param[2] = abs(fit.param[2])
    return fit, model
end
function fit_corr_2exp(c,w,ncut;sign=+1)
    fit_1exp = fit_corr_1exp(c,w,ncut;sign)[1]
    m0, f0 = fit_1exp.param[1], fit_1exp.param[2]
    T = length(c)
    t = collect(1:T)
    cut = _aply_cut(ncut,T)
    if ndims(w)==2
        weight = w[cut,cut]
        weight = (weight + weight')/2
    else
        weight = w[cut]
    end
    @. model(t,p) = model_1st(t-1,T,p[1],p[2];sign) + model_2nd(t-1,T,p[3],p[4];sign)
    p0 = [m0,f0,1.1*m0,f0*f0*m0]
    fit = curve_fit(model,t[cut],c[cut],weight,p0)
    if fit.param[1] > fit.param[3] || fit.param[2] < 10^-12
        fit.param[1], fit.param[3] = fit.param[3], fit.param[1]
        fit.param[2], fit.param[4] = fit.param[4], fit.param[2]
    end
    return fit, model
end
function fit_corr(c,cov,ncut;constant=false,nexp2=false,sign=+1)
    if constant
        return fit_corr_const(c,cov,ncut;sign)
    end
    if ndims(cov)==2
        w  = inv(cov)
        w = (w + w')/2
    else
        w = inv.(cov)
    end
    return nexp2 ? fit_corr_2exp(c,w,ncut;sign) : fit_corr_1exp(c,w,ncut;sign)
end
function fit_corr_bars(c,cov,ncut;nexp2=false,constant=false,sign=+1)
    if ndims(cov)==2
        Δc = sqrt.(diag(cov))
    else
        Δc = sqrt.(cov)
    end
    if constant
        fit1, fit2, fit3, model = fit_corr_bars_const(c,cov,ncut;sign)
    else
        fit1, model = fit_corr(c,cov,ncut;nexp2,sign)
        fit2, model = fit_corr(c+Δc,cov,ncut;nexp2,sign)
        fit3, model = fit_corr(c-Δc,cov,ncut;nexp2,sign)
    end
    return fit1, fit2, fit3, model
end
###############################################################
# fit additional constant
# automatically fits a constant present in the correlator
# that can appear in multiple circumstances
# i) residual constant in σ meson 2) topological effect in η'
# meson 3) around-the-world effects in 2-meson-states
###############################################################
function fit_corr_const(c,cov,ncut;sign=+1)
    T = length(c)
    t = collect(1:T)
    cut = _aply_cut(ncut,T)
    if ndims(cov)==2
        w = inv(cov)
        w = (w + w')/2
        weight = w[cut,cut]
        weight = (weight + weight')/2
    else
        w = inv.(cov)
        weight = w[cut]
    end

    # give bounds, otherwise LsqFit likes to fit to only the constant
    lower = [0.0,-Inf64, -2maximum(abs.(c))]
    upper = [2.0, Inf64, +2maximum(abs.(c))]
    initial = [0.1,1.0,0.0]

    @. model(t,p) = model_2nd(t-1,T,p[1],p[2];sign) + p[3]
    fit = curve_fit(model,t[cut],c[cut],weight,initial;lower,upper)
    fit.param[2] = abs(fit.param[2])
    return fit, model
end
function fit_corr_bars_const(c,cov,ncut;sign=+1)
    if ndims(cov)==2
        Δc = sqrt.(diag(cov))
    else
        Δc = sqrt.(cov)
    end
    fit1, model = fit_corr_const(c,cov,ncut;sign)
    fit2, model = fit_corr_const(c+Δc,cov,ncut;sign)
    fit3, model = fit_corr_const(c-Δc,cov,ncut;sign)
    return fit1, fit2, fit3, model
end
