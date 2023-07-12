module HiRepAnalysis

import LsqFit: curve_fit, stderror, confidence_interval, coef
using Plots
using Statistics
using LinearAlgebra
using Parsers
using Distributions
using ProgressMeter
using NaNStatistics
using HDF5
using Roots
gr(legendfontsize=12,tickfontsize=12,ms=7,msw=2)

include("average.jl")
export average_plaquette, average_correlator
include("fitcorr.jl")
export fitmass, fit_corr, fit_corr_bars
include("mass.jl")
export effectivemass, effectivemass_err, effectivemass_cosh, effectivemass_cosh_err, implicit_meff
include("parse.jl")
export HMC_accept, plaquettes, correlators, latticesize, loweig, quarkmasses, gaugegroup, couplingβ, fileinfo, confignames
export dilution, ncolors, nconfigs, nhits, parse_disconnected
include("disconnected.jl")
export flatten_disc,  hit_time_average_disconnected, disconnected_eta, disconnected_pi0, disconnected_eta_MC, disconnected_pi0_MC
include("mesonanalysis.jl")
export meson_mass_decay, meson_mass_decay_select, meson_mass_decay_bootstrap, meson_mass_decay_jackknife, correlator_deriv
export apply_jackknife, singlet_jackknife, apply_bootstrap, singlet_bootstrap, groundstate_correlator
include("errorstring.jl")
export errorstring
include("autocor.jl")
export autocorrelation_time, autocorrelation
include("pcac.jl")
export awi_fit, awi_correlator, awi_mass, nondeg_pcac_mass
include("histogramfit.jl")
export decay_mass_histogram

end # module
