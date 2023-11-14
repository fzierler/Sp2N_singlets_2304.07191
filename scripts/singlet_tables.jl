using DelimitedFiles
using LaTeXStrings
function write_tex_table(name,data;insert_hline=2,labels="")
    io = open(name,"w")
    rows, cols = size(data)
    # I have hardcoded centering of the table's contents
    table_layout = repeat("|c",cols)*"|"
    header = """\\begin{tabular}{$table_layout}
    \t\\hline
    \t$labels\\\\
    \t\\hline\\hline
    """
    write(io,header)
    for i in 1:rows
        # if the gauge coupling changes insert two hlines for formatting
        if i > 1 && data[i,insert_hline] != data[i-1,insert_hline]
            write(io,"\t\\hline\\hline\n")
        end
        # we could use some padding here for nicer formatting but for now
        # I only insert a '&' to make a minimal tex-compliant table
        write(io,"\t"*string(data[i,1])*"&")
        for j in 2:cols-1
            write(io,string(data[i,j])*"&")
        end
        write(io,string(data[i,cols])*"\\\\\n")
    end
    write(io,"\t\\hline\\hline\n")
    write(io,"\\end{tabular}")
    close(io)
end
function generate_singlet_tables()
    # obtain gauge group for every ensemble
    group = readdlm("input/parameters/param_deriv.csv",';',skipstart=1)[:,10]
    # read in the tabulated results

    # MAIN PAPER: PARAMETERS
    # Main results for the singlet mesons with degenerate fermions, TABLE I
    data_v0 = readdlm("output/data_human_readable/data_deriv_HR.csv",';',skipstart=1)
    param_v0 = readdlm("input/parameters/param_deriv.csv",';',skipstart=1)

    p0 = param_v0[:,11]
    p1 = data_v0[:,1:11]
    p4 = data_v0[:,17]

    data = hcat(p0,group,p1,p4)
    data = replace(data,NaN => "-")
    data = replace(data,-1 => "-")
    #writedlm("output/tables/tab_ensembles_degenerate.csv",data,"&")
    labels=L"Ensemble & group & $\beta$ & $m_0$ & $L$ & $T$ & $n_\text{conf}$ & $n_\text{src}$ & $I_{\eta'}$ & $I_\pi$ & $I_\sigma$ & $I_{\sigma^{\rm conn.}}$ & $I_\rho$ &  $\langle P \rangle$"
    isdir("output/tables") || mkdir("output/tables")
    write_tex_table("output/tables/tab_ensembles_degenerate.tex",data;insert_hline=3,labels)

    # NON-DEG PARAMETERS
    # Main results for the singlet mesons with degenerate fermions, TABLE II
    data_v0 = readdlm("output/data_human_readable/data_non_deg_HR.csv",';',skipstart=1)
    param_v0 = readdlm("input/parameters/param_non_deg.csv",';',skipstart=1)

    data = hcat(param_v0[:,10],data_v0[:,1:14])
    data = replace(data,-1 => "-")
    data = replace(data,NaN => "-")
    #writedlm("output/tables/tab_ensembles_non-degenerate.csv",data,"&")
    labels=L"Ensemble & $\beta$ & $m_0^1$ & $m_0^2$ & $L$ & $T$ & $n_\text{conf}$ & $n_\text{src}$ & $I_{\eta'}$ & $I_{\pi^0}$ & $I_{\pi^\pm}$ & $I_\sigma$ & $I_{\sigma^{\rm conn.}}$ & $I_\rho$ &  $\langle P \rangle$"
    write_tex_table("output/tables/tab_ensembles_non-degenerate.tex",data[1:end-1,:];labels)

    # MAIN PAPER: RESULTS
    # Main results for the singlet mesons with degenerate fermions, TABLE III
    data_v0 = readdlm("output/data_human_readable/data_deriv_HR.csv",';',skipstart=1)

    ens = data_v0[:,1:4]
    p1 = data_v0[:,19]
    p2 = data_v0[:,18]
    p3 = data_v0[:,12:13]
    p4 = data_v0[:,15:16]

    data = hcat(group,ens,p1,p2,p3,p4)
    data = replace(data,NaN => "-")
    data = replace(data,-1 => "-")
    #writedlm("output/tables/tab_masses_degenerate.csv",data,"&")
    labels=L"& $\beta$ & $m_0$ & $L$ & $T$ & $m_\pi L$ & $m_\pi / m_\rho$ & $m_\pi$ & $m_\rho$ & $m_{\eta'}$ & $m_\sigma$"
    write_tex_table("output/tables/tab_masses_degenerate.tex",data;labels)

    # NON-DEG RESULTS
    # Main results for the singlet mesons with non-degenerate fermions, TABLE IV
    data_v0 = readdlm("output/data_human_readable/data_non_deg_HR.csv",';',skipstart=1)
    p1 = data_v0[:,1:5]
    p2 = data_v0[:,15:19]
    p3 = data_v0[:,21]
    p4 = data_v0[:,20]
    data = hcat(p1,p2,p3,p4)

    data = replace(data,NaN => "-")
    data = replace(data,-1 => "-")
    #writedlm("output/tables/tab_masses_non-degenerate.csv",data,"&")
    labels=L"$\beta$ & $m_0^{(1)}$ & $m_0^{(2)}$ & $L$ & $T$ & $m_{\pi^0}$ & $m_{\pi^0_c}$ & $m_{\pi^\pm}$ & $m_{\rho^0}$ & $m_{\rho^\pm}$ & $m_{\eta'}$ & $m_{\sigma}$"
    write_tex_table("output/tables/tab_masses_non-degenerate.tex",data[1:end-1,:];labels)

    # APPENDIX: ETA MESON
    # Comparison of different analysis methods, TABLE V
    data_v0 = readdlm("output/data_human_readable/data_deriv_HR.csv",';',skipstart=1)
    data_v1 = readdlm("output/data_human_readable/data_vsub_HR.csv",';',skipstart=1)
    data_v2 = readdlm("output/data_human_readable/data_nothing_HR.csv",';',skipstart=1)
    data_v3 = readdlm("output/data_human_readable/data_const_HR.csv",';',skipstart=1)
    data_v4 = readdlm("output/data_human_readable/data_only_deriv_HR.csv",';',skipstart=1)

    ens = data_v0[:,1:4]
    mη0 = data_v0[:,15]
    mη1 = data_v1[:,15]
    mη2 = data_v2[:,15]
    mη3 = data_v3[:,15]
    mη4 = data_v4[:,15]

    data = hcat(group,ens,mη0,mη1,mη2,mη3,mη4)
    data = replace(data,NaN => "-")
    #writedlm("output/tables/tab_eta_different_techniques.csv",data,"&")
    labels=L"& $\beta$ & $m_0$ & $L$ & $T$ & $m_{\eta'}$ & $m_{\eta'}^{\rm (i)}$ & $m_{\eta'}^{\rm (ii)}$ & $m_{\eta'}^{\rm (iii)}$ & $m_{\eta'}^{\rm (iv)}$"
    write_tex_table("output/tables/tab_eta_different_techniques.tex",data;labels)

    # APPENDIX: SIGMA MESON
    # Comparison of different analysis methods, TABLE VI
    data_v4 = readdlm("output/data_human_readable/data_deriv_HR.csv",';',skipstart=1)
    data_v5 = readdlm("output/data_human_readable/data_both_HR.csv",';',skipstart=1)
    data_v6 = readdlm("output/data_human_readable/data_vsub_HR.csv",';',skipstart=1)
    data_v7 = readdlm("output/data_human_readable/data_both_nogs_HR.csv",';',skipstart=1)

    ens = data_v0[:,1:4]
    mσ0 = data_v4[:,16]
    mσ1 = data_v5[:,16]
    mσ2 = data_v6[:,16]
    mσ3 = data_v7[:,16]

    data = hcat(group,ens,mσ0,mσ1,mσ2,mσ3)
    data = replace(data,NaN => "-")
    #writedlm("output/tables/tab_sigma_different_techniques.csv",data,"&")
    labels=L"& $\beta$ & $m_0$ & $L$  & $T$  & $m_{\sigma}$ & $m_{\sigma}^{(i)}$ & $m_{\sigma}^{\rm (ii)}$ & $m_{\sigma}^{\rm (iii)}$"
    write_tex_table("output/tables/tab_sigma_different_techniques.tex",data;labels)
end
generate_singlet_tables()
