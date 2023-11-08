using HiRepAnalysis
using LinearAlgebra
using DelimitedFiles

parse_fit_intervall(s) = isa(s,Integer) ? s : Tuple(parse.(Int,split(s[2:end-1],',')))
function read_prm_nondeg(file,i)
    prm = readdlm(file,';')
    odir = prm[i,1]
    hits = prm[i,2]
    fitη, fitπ0, fitπc, fitσ, fita0, fitρ = parse_fit_intervall.(prm[i,3:8])
    name = prm[i,10]
    return odir, hits, fitη, fitπ0, fitπc, fitσ, fita0, fitρ, name
end
function fermion_masses_from_filename(file)
    p1 = last(findfirst("m1",file))
    p2 = first(findfirst("m2",file))
    p3 = last(findfirst("m2",file))
    # create array of masses for matching of output
    m  = [file[p1+1:p2-1],file[p3+1:end]]
    return m
end

prmfile = "input/parameters/param_non_deg.csv"
@showprogress for i in 2:countlines(prmfile)
    odir, hits, fitη, fitπ0, fitπc, fitσ, fita0, fitρ, ensemble_name = read_prm_nondeg(prmfile,i)
    m = fermion_masses_from_filename(odir)

    type    = "DISCON_SEMWALL SINGLET"
    typeM   = "DEFAULT_SEMWALL TRIPLET"
    typeU   = "SEMWALL_U TRIPLET"
    typeD   = "SEMWALL_D TRIPLET"
    typeUD  = "SEMWALL_UD TRIPLET"
    fileM   = joinpath(odir,"out_spectrum_with_pcac.hdf5")
    fileS   = joinpath(odir,"out_spectrum_discon.hdf5")
    hdf5groupM = splitext(fileM)[1]
    hdf5groupS = splitext(fileS)[1]

    T = latticesize(hdf5file,hdf5groupM)[1]
    L = latticesize(hdf5file,hdf5groupM)[2]
    rescale = (L^3)^2 /L^3
    writefile = true

    # obtain disconnected part
    if m[1] != m[2]
        C_dis_MC_sigma = disconnected_eta_MC(hdf5file,hdf5groupS,type,hits,m[1],m[2];rescale=rescale,key="id_disc_re",vsub=false)
        # obtain connected part
        C_con_MC = correlators(hdf5file,hdf5groupM,typeU,typeD,"g5")
        C_con_MC_sigma = correlators(hdf5file,hdf5groupM,typeU,typeD,"id")
        HiRepAnalysis._rescale_corrs!(C_con_MC, L)
        # obtain non-singlet masses
        kws = (error=:jack,nexp2=false)
        mπc, Δmπc = meson_mass_decay(hdf5file,hdf5groupM,"g5",typeU,typeD;ncut=fitπc,kws...)[1:2]
        mρc, Δmρc = meson_mass_decay(hdf5file,hdf5groupM,"g1",typeU,typeD;ncut=fitρ,kws...)[1:2]
        mπf, Δmπf = meson_mass_decay(hdf5file,hdf5groupM,"g5",typeUD;ncut=fitπc,kws...)[1:2]
        mρf, Δmρf = meson_mass_decay(hdf5file,hdf5groupM,"g1",typeUD;ncut=fitρ,kws...)[1:2]
        #obtain PCAC/AWI mass
        kwsAWI = (therm=0,step=1,autocor=false,ncut=fitπc[1])
        mAWI, ΔmAWI = awi_mass(hdf5file,hdf5groupM,typeUD;kwsAWI...)
        # Read in correct values of the variational analysis
        eta_pi_filename = "output/data/non_degenerate_data_eta_pi.csv"
        eta_pi_data = readdlm(eta_pi_filename,';')
        # compare filename and find correct row:
        filenames = eta_pi_data[:,1]
        ind = findfirst(contains(basename(odir)),filenames)
        if fitπ0 != -1
            mπ0, Δmπ0 = eta_pi_data[ind,7], eta_pi_data[ind,8]
        else
            mπ0, Δmπ0 = NaN,NaN
        end
        if fitη != -1
            mη,  Δmη  = eta_pi_data[ind,9], eta_pi_data[ind,10]
        else
            mη,  Δmη = NaN,NaN
        end
    else
        C_dis_MC_eta = disconnected_eta_MC(hdf5file,hdf5groupS,type,hits;rescale=rescale,key="g5_disc_re",vsub=false)
        C_dis_MC_pi0 = zero(C_dis_MC_eta)
        C_dis_MC_sigma = disconnected_eta_MC(hdf5file,hdf5groupS,type,hits;rescale=rescale,key="id_disc_re",vsub=false)
        # obtain connected part
        C_con_MC = correlators(hdf5file,hdf5groupM,typeM,"g5")
        C_con_MC_sigma = correlators(hdf5file,hdf5groupM,typeM,"id")
        HiRepAnalysis._rescale_corrs!(C_con_MC, L)
        #obtain PCAC/AWI mass
        kwsAWI = (therm=0,step=1,autocor=false,ncut=fitπc[1])
        mAWI, ΔmAWI = awi_mass(hdf5file,hdf5groupM,typeM;kwsAWI...)
        # obtain non-singlet masses
        kws = (error=:jack,nexp2=false)
        mπf, Δmπf = meson_mass_decay(hdf5file,hdf5groupM,"g5",typeM;ncut=fitπc,kws...)[1:2]
        mρf, Δmρf = meson_mass_decay(hdf5file,hdf5groupM,"g1",typeM;ncut=fitρ,kws...)[1:2]
        mπc, Δmπc = mπf, Δmπf
        mρc, Δmρc = mρf, Δmρf
        # perform analysis only for degenerate masses
        if size(C_dis_MC_eta)[1] != size(C_con_MC)[2]
            @warn "mismatch between connected  $(size(C_con_MC)[2]) and disconnected $(size(C_dis_MC_eta)[1])"
        end
        if fitπ0 != -1
            mπ0, Δmπ0, Cdπ0 ,ΔCdπ0, Cπ0, ΔCπ0, meffπ0, Δmeffπ0, cπ0, Δcπ0 = singlet_jackknife(C_con_MC,C_dis_MC_pi0,fitπc,fitπ0;deriv=true,gs_sub=true,sigma=false,constant=false)
        else
            mπ0, Δmπ0 = NaN,NaN
        end
        if fitη != -1
            mη, Δmη, Cdη, ΔCdη, Cη, ΔCη, meffη, Δmeffη, cη, Δcη = singlet_jackknife(C_con_MC,C_dis_MC_eta,fitπc,fitη;deriv=true,gs_sub=true,sigma=false,constant=false)
        else
            mη,  Δmη = NaN,NaN
        end
    end
    if fitσ != -1
        mσ, Δmσ, Cdσ, ΔCdσ, Cσ, ΔCσ, meffσ, Δmeffσ, cσ, Δcσ = singlet_jackknife(C_con_MC_sigma,C_dis_MC_sigma,fita0,fitσ;deriv=true,gs_sub=true,sigma=true,constant=false)
    else
        mσ, Δmσ = NaN,NaN
    end

    G = gaugegroup(hdf5file,hdf5groupM)
    β = couplingβ(hdf5file,hdf5groupM)
    maxconf = nconfigs(hdf5file,hdf5groupM)
    name="$(T)x$(L)^3 $G, m=$(first(m)),$(last(m))"

    P, ΔP = average_plaquette(hdf5file,hdf5groupM)

    sπc = errorstring(mπc,Δmπc)
    sρc = errorstring(mρc,Δmρc)
    sπf = errorstring(mπf,Δmπf)
    sρf = errorstring(mρf,Δmρf)
    sη  = errorstring(mη,Δmη)
    sσ  = errorstring(mσ,Δmσ)
    sπ0 = errorstring(mπ0,Δmπ0)
    # calculate ratios
    Δratio(x,y,Δx,Δy) = sqrt((Δx / y)^2 + (Δy*x / y^2)^2)
    # more errorstrings
    sP  = errorstring(P,ΔP)
    sAWI = errorstring(mAWI,ΔmAWI)

    if writefile
        isdir("output/data/") || mkdir("output/data/")
        isdir("output/data_human_readable/") || mkdir("output/data_human_readable/")
        ioPR = open("output/data/data_non_deg.csv","a+")
        iszero(position(ioPR)) && write(ioPR,"β;m1;m2;L;T;Nconf;hits;fitη;fitπ0;fitπc;fitσ;fita0;fitρ;P;ΔP;mπ0;Δmπ0;mπc;Δmπc;mπf;Δmπf;mρ0;Δmρ0;mρf;Δmρf;mσ;Δmσ;mη;Δmη;mAWI;ΔmAWI\n")
        write(ioPR,"$β;$(m[1]);$(m[2]);$L;$T;$maxconf;$hits;$fitη;$fitπ0;$fitπc;$fitσ;$fita0;$fitρ;$P;$ΔP;$mπ0;$Δmπ0;$mπc;$Δmπc;$mπf;$Δmπf;$mρc;$Δmρc;$mρf;$Δmρf;$mσ;$Δmσ;$mη;$Δmη;$mAWI;$ΔmAWI\n")
        close(ioPR)
        ioPR = open("output/data_human_readable/data_non_deg_HR.csv","a+")
        iszero(position(ioPR)) && write(ioPR,"β;m1;m2;L;T;Nconf;hits;fitη;fitπ0;fitπc;fitσ;fita0;fitρ;P;mπ0;mπc;mπf;mρ0;mρf;mσ;mη;mAWI\n")
        write(ioPR,"$β;$(m[1]);$(m[2]);$L;$T;$maxconf;$hits;$fitη;$fitπ0;$fitπc;$fitσ;$fita0;$fitρ;$sP;$sπ0;$sπc;$sπf;$sρc;$sρf;$sσ;$sη;$sAWI\n")
        close(ioPR)
    end
end
