using DelimitedFiles
using Plots
using LaTeXStrings
using Plots.PlotMeasures
using ColorSchemes
pgfplotsx(ms=6,legendfontsize=16,titlefontsize=20,labelfontsize=20,tickfontsize=20,framestyle=:box)

file ="input/su3_literature.csv"
data = readdlm(file,'\t')[2:end,:]

mπ = Float64.(data[:,3])
Δmπ= Float64.(data[:,4])
mρ = Float64.(data[:,5])
Δmρ= Float64.(data[:,6])
mη = Float64.(data[:,7])
Δmη= Float64.(data[:,8])
r  = Float64.(data[:,9])
Δr = Float64.(data[:,10])
collab = String.(data[:,1])

Δratio(x,y,Δx,Δy) = sqrt((Δx / y)^2 + (Δy*x / y^2)^2)
mρπ = @. mρ / mπ
mηπ = @. mη / mπ
mηρ = @. mη / mρ
Δmρπ = Δratio.(mρ,mπ,Δmρ,Δmπ)
Δmηπ = Δratio.(mη,mπ,Δmη,Δmπ)
Δmηρ = Δratio.(mη,mρ,Δmη,Δmρ)

collaboration_indices  = [1:6, 7:7, 8:9, 10:11, 12:23, 24:26, 27:27, 29:32, 33:35]
colors   = [:white, :grey]
colors   = palette(:Dark2_3)
markers  = deleteat!(Plots.supported_markers(),(1,2,4,7,9,13,14,16,17,18))
append!(markers,[:circle])
mss = [6,8,7,8,7,6,7,7,6] # markersiszes

plt1 = plot(title=L"$\rm{SU}(3)$ with $N_f=2$",xlabel=L"m_\pi / m_\rho",ylabel=L" m_{\rm meson} / m_\pi",legend=:topright)
# legend for eta' and rho
scatter!(plt1,[0],[0],markershape=:circle,mc=colors[1],label=L"\rho")
scatter!(plt1,[0],[0],markershape=:circle,mc=colors[2],label=L"\eta'")
# plot data
for (n,i) in enumerate(collaboration_indices)
    scatter!(plt1,r[i],xerr=Δr[i],mηπ[i],yerr=Δmηπ[i],color=colors[2],ms=mss[n],marker=markers[n],label="")
    scatter!(plt1,r[i],xerr=Δr[i],mρπ[i],yerr=Δmρπ[i],color=colors[1],ms=mss[n],marker=markers[n],label="")
end
# now create legends
for (n,i) in enumerate(collaboration_indices)
    scatter!(plt1,[0],[0],color=:white,ms=mss[n],marker=markers[n],label=first(collab[i]))
end
plot!(plt1,xlims=(0.12,0.9),ylims=(0.8,7))
plot!(plt1,legendfontsize=14)

# plot data
colors = palette(:Set1_9)
plt2 = plot(title=L"$\rm{SU}(3)$ with $N_f=2$",xlabel=L"m_\pi / m_\rho",ylabel=L" m_{\eta'} / m_\rho ",legend=:outerright)
for (n,i) in enumerate(collaboration_indices)
    scatter!(plt2,r[i],xerr=Δr[i],mηρ[i],yerr=Δmηρ[i],color=colors[n],ms=mss[n],marker=markers[n],label="")
end
# now create legends
for (n,i) in enumerate(collaboration_indices)
    scatter!(plt2,[0],[0],color=colors[n], ms=mss[n],marker=markers[n],label=first(collab[i]))
end
plot!(plt2,xlims=(0.12,0.9),ylims=(0.55,1.35))


colors = palette(:Dark2_3)
# plots for comparing gauge groups
plt3 = scatter(r,mηρ,yerr=Δmηρ,xerr=Δr,title="",legend=:outerright,color=colors[1],label=L"\rm{SU(3)}",ms=8,markershape=:diamond)
plt4 = scatter(r,mηπ,yerr=Δmηπ,xerr=Δr,title="",legend=:topright,color=colors[1],label=L"\rm{SU(3)}",ms=8,markershape=:diamond)

# Now add our SU(2) and Sp(4) data
data = readdlm("output/data/data_deriv.csv",';',skipstart=1)
β,mq,L,T = ( data[:,i] for i in 1:4 )
mπ,Δmπ,mρ,Δmρ,ma0,Δma0,mη,Δmη,mσ,Δmσ = ( data[:,i] for i in 12:21 )
# get mass ratios
x   = @. mπ / mρ
mηπ = @. mη / mπ
mηρ = @. mη / mρ
Δx   = Δratio.(mπ,mρ,Δmπ,Δmρ)
Δmηπ = Δratio.(mη,mπ,Δmη,Δmπ)
Δmηρ = Δratio.(mη,mρ,Δmη,Δmρ)
# get indices of comparable sets of ensembles_degenerate
iβ1 = findall(x -> isapprox(6.9,x),β)
iβ2 = findall(x -> isapprox(7.2,x),β)
iSU2   = findall(x -> isapprox(2.0,x),β)

scatter!(plt3,x[iSU2],mηρ[iSU2],yerr=Δmηρ[iSU2],xerr=Δx[iSU2],color=colors[2],label=L"\rm{SU(2)}(\beta=2.0)",ms=8,markershape=:pentagon)
scatter!(plt4,x[iSU2],mηπ[iSU2],yerr=Δmηπ[iSU2],xerr=Δx[iSU2],color=colors[2],label=L"\rm{SU(2)}(\beta=2.0)",ms=8,markershape=:pentagon)
scatter!(plt3,x[iβ1],mηρ[iβ1],yerr=Δmηρ[iβ1],xerr=Δx[iβ1],color=colors[3],label=L"\rm{Sp(4)}(\beta=6.9)",ms=8,markershape=:rect)
scatter!(plt4,x[iβ1],mηπ[iβ1],yerr=Δmηπ[iβ1],xerr=Δx[iβ1],color=colors[3],label=L"\rm{Sp(4)}(\beta=6.9)",ms=8,markershape=:rect)
scatter!(plt3,x[iβ2],mηρ[iβ2],yerr=Δmηρ[iβ2],xerr=Δx[iβ2],color=colors[3],label=L"\rm{Sp(4)}(\beta=7.2)",ms=8,markershape=:circle)
scatter!(plt4,x[iβ2],mηπ[iβ2],yerr=Δmηπ[iβ2],xerr=Δx[iβ2],color=colors[3],label=L"\rm{Sp(4)}(\beta=7.2)",ms=8,markershape=:circle)
plot!(plt3,xlims=(0.63,0.9),ylims=(0.55,1.42),xlabel=L"m_\pi / m_\rho",ylabel=L" m_{\eta'} / m_\rho ")
plot!(plt4,xlims=(0.63,0.9),ylims=(0.98,1.50),xlabel=L"m_\pi / m_\rho",ylabel=L" m_{\eta'} / m_\pi ")

plt_combined = plot(plt1,plt2,plt4,plt3,layout=(2,2),size=(1.6*700,1.6*500),right_margin = 10mm, bottom_margin = 10mm)
savefig(plt_combined,"output/figures/SU3_vs_Sp4.pdf")
plt_combined