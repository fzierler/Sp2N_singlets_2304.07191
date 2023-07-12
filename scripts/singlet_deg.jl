using HiRepAnalysis
using LinearAlgebra
using DelimitedFiles

parse_fit_intervall(s) = isa(s,Integer) ? s : Tuple(parse.(Int,split(s[2:end-1],',')))
function read_prm(file,i)
    prm = readdlm(file,';')
    odir = prm[i,1]
    hits = prm[i,2]
    fitη,fitπ,fitσ,fita0,fitρ = parse_fit_intervall.(prm[i,3:7])
    fileM, fileS = joinpath(odir,"out_spectrum.hdf5"), joinpath(odir,prm[i,9])
    group,name,vsub,deriv,gs_sub,constant = prm[i,10:15]
    return odir, hits, fitη, fitπ, fitσ, fita0, fitρ ,fileM, fileS, group, name, vsub, deriv, gs_sub, constant
end
function singlet_analysis(hdf5file,hdf5groupM,hdf5groupS,hits,vsub,key,cut,fitint,deriv,gs_sub,sigma,constant;maxconf=typemax(Int))
    type  = "DISCON_SEMWALL SINGLET"
    typeM = "DEFAULT_SEMWALL TRIPLET"
    T, L = latticesize(hdf5file,hdf5groupM)[1:2]
    rescale = (L^3)^2 /L^3
    # disconnected contributions from code in /Spectrum/
    # perform correct normalization
    C_dis_MC = disconnected_eta_MC(hdf5file,hdf5groupS,type,hits;maxhits=hits,rescale,vsub,key="$(key)_disc_re")
    C_con_MC = correlators(hdf5file,hdf5groupM,typeM,key)
    HiRepAnalysis._rescale_corrs!(C_con_MC, L)
    # restrict analysis to first N configurations
    if size(C_dis_MC)[1] != size(C_con_MC)[2]
        @warn "mismatch between connected  $(size(C_con_MC)[2]) and disconnected $(size(C_dis_MC)[1])"
        N = min(size(C_dis_MC)[1],size(C_con_MC)[2],maxconf)
        C_dis_MC = C_dis_MC[1:N,:]
        C_con_MC = C_con_MC[:,1:N]
    end
    mη, Δmη, Cd,ΔCd, C, ΔC, meffη, Δmeffη, cη, Δcη = singlet_jackknife(C_con_MC,C_dis_MC,cut,fitint;deriv,gs_sub,sigma,constant)
    return mη, Δmη, meffη, Δmeffη, cη, Δcη, C, ΔC, Cd, ΔCd
end
function write_singlet_files(filedat,fileM,fitη,fitπ,fitσ,fita0,fitρ,mπ,Δmπ,mρ,Δmρ,ma0,Δma0,mη,Δmη,mσ,Δmσ,hits,P,ΔP,T,L,Nconf,m,G,β)
    # calculate ratios
    Δratio(x,y,Δx,Δy) = sqrt((Δx / y)^2 + (Δy*x / y^2)^2)
    πρ  = mπ / mρ
    πL  = L*mπ
    Δπρ = Δratio(mπ,mρ,Δmπ,Δmρ)
    ΔπL = L*Δmπ
    # and human readable strings
    sπρ = errorstring(πρ,Δπρ)
    sπ  = errorstring(mπ,Δmπ)
    sρ  = errorstring(mρ,Δmρ)
    sη  = errorstring(mη,Δmη)
    sσ  = errorstring(mσ,Δmσ)
    sa0 = errorstring(ma0,Δma0)
    sπL = errorstring(πL,ΔπL)
    sP  = errorstring(P,ΔP)
    # set up csv's for
    isdir("output/data/") || mkdir("output/data/")
    isdir("output/data_human_readable/") || mkdir("output/data_human_readable/")
    ioPR = open("output/data/"*filedat*".csv","a+")
    iszero(position(ioPR)) && write(ioPR,"β;m;L;T;Nconf;Nhits;fitη;fitπ;fit_sigma;fit_a0;fit_rho;mπ;Δmπ;mρ;Δmρ;ma0;Δma0;mη;Δmη;mσ;Δmσ;P;ΔP;mπ/mρ;Δ(mπ/mρ);mπL;ΔmπL\n")
    write(ioPR,"$β;$(m[1]);$L;$T;$Nconf;$hits;$fitη;$fitπ;$fitσ;$fita0;$fitρ;$mπ;$Δmπ;$mρ;$Δmρ;$ma0;$Δma0;$mη;$Δmη;$mσ;$Δmσ;$P;$ΔP;$πρ;$Δπρ;$πL;$ΔπL\n")
    close(ioPR)
    ioPR = open("output/data_human_readable/"*filedat*"_HR.csv","a+")
    iszero(position(ioPR)) && write(ioPR,"β;m;L;T;Nconf;hits;fitη;fitπ;fit_sigma;fit_a0;fit_rho;mπ;mρ;ma0;mη;mσ;P;mπ/mρ;mπL\n")
    write(ioPR,"$β;$(m[1]);$L;$T;$Nconf;$hits;$fitη;$fitπ;$fitσ;$fita0;$fitρ;$sπ;$sρ;$sa0;$sη;$sσ;$sP;$sπρ;$sπL\n")
    close(ioPR)
end

files = readdir("input/parameters/";join=true)
filter!(!contains("flow"),files)
filter!(!contains("non_deg"),files)
N_ensembles = (countlines(files[1])-1)*length(files)
p = Progress(N_ensembles)
for prmfile in files
    for i in 2:countlines(prmfile)
        id = replace(first(splitext(basename(prmfile))),"param_" => "")
        odir, hits, fitη, fitπ, fitσ, fita0, fitρ, fileM, fileS, group, name, vsub, deriv, gs_sub, constant = read_prm(prmfile,i)

        minplateau = 4
        hdf5groupM = splitext(fileM)[1]
        hdf5groupS = splitext(fileS)[1]
        Nconf = nconfigs(hdf5file,hdf5groupM)
        m = quarkmasses(hdf5file,hdf5groupM)
        G = gaugegroup(hdf5file,hdf5groupM)
        β = couplingβ(hdf5file,hdf5groupM)
        T, L = latticesize(hdf5file,hdf5groupM)[1:2]

        kws = (therm=0,step=1,autocor=true)
        typeM = "DEFAULT_SEMWALL TRIPLET"
        mπ, Δmπ = meson_mass_decay_select(hdf5file,hdf5groupM,"g5",typeM;error=:jack,nexp2=false,ncut=fitπ  ,kws...)[1:2]
        mρ, Δmρ = meson_mass_decay_select(hdf5file,hdf5groupM,"g1",typeM;error=:jack,nexp2=false,ncut=fitρ  ,kws...)[1:2]
        if fita0 != -1
            ma, Δma = meson_mass_decay_select(hdf5file,hdf5groupM,"id",typeM;error=:jack,nexp2=false,ncut=fita0 ,kws...)[1:2]
        else
            ma, Δma = NaN, NaN
        end

        T, L = latticesize(hdf5file,hdf5groupM)[1:2]
        rescale = (L^3)^2 /L^3
        type = "DISCON_SEMWALL SINGLET"
        P, ΔP = average_plaquette(hdf5file,hdf5groupM)

        if fitη != -1
            key, sigma = "g5", false
            mη, Δmη, meffη, Δmeffη, cη, Δcη, C, ΔC, Cd, ΔCd = singlet_analysis(hdf5file,hdf5groupM,hdf5groupS,hits,vsub,key,fitπ,fitη,deriv,gs_sub,sigma,constant)
        end
        if fitσ != -1
            if fita0 == -1
                gs_sub = false
            end
            key, sigma, constant = "id", true, false
            mσ, Δmσ, meffσ, Δmeffσ, cσ, Δcσ, C, ΔC, Cd, ΔCd = singlet_analysis(hdf5file,hdf5groupM,hdf5groupS,hits,vsub,key,fita0,fitσ,deriv,gs_sub,sigma,constant)
        end

        if (fitη == -1) || (fitη[2] - fitη[1] < minplateau - 1)
            mη, Δmη = NaN, NaN
        end
        if (fitσ == -1) || (fitσ[2] - fitσ[1] < minplateau - 1)
            mσ, Δmσ = NaN, NaN
        end
        write_singlet_files("data_$id",fileM,fitη,fitπ,fitσ,fita0,fitρ,mπ,Δmπ,mρ,Δmρ,ma,Δma,mη,Δmη,mσ,Δmσ,hits,P,ΔP,T,L,Nconf,m,G,β)
        next!(p)
    end
end
