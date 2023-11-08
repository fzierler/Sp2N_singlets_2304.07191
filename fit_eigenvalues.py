import gvar as gv
import numpy as np
import corrfitter as cf
import h5py as h5
import csv
import os
from ast import literal_eval

def make_models(T,tmin,tmax):
    tdata = range(T)
    tfit = tdata[tmin:tmax]
    models = [cf.Corr2(datatag='Gaa', tdata=tdata, tfit=tfit, a='a', b='a', dE='dE',tp=-T)]
    return models

def make_prior(N):
    prior = gv.BufferDict()
    prior['log(a)']  = gv.log(gv.gvar(N * ['1(1)']))
    prior['log(dE)'] = gv.log(gv.gvar(N * ['1(1)']))
    return prior

def main_fitting(data,N,T,tmin,tmax):
    models = make_models(T,tmin,tmax)
    prior = make_prior(N)
    fitter = cf.CorrFitter(models=models)
    fit = fitter.lsqfit(data=data, prior=prior)
    return fit

def main(outfileHR, outfile ,h5file,Nexp, tmin_pi, tmax_pi, tmin_eta, tmax_eta ):
    fname = h5file
    f = h5.File(fname,'r')

    eigvals_deriv = f['eigvals_deriv']
    Delta_eigvals_deriv = f['Δeigvals_deriv']
    T = f['lattice'][0]
    
    eig1 = dict(Gaa=gv.gvar(eigvals_deriv[:,0],Delta_eigvals_deriv[:,0]))
    eig2 = dict(Gaa=gv.gvar(eigvals_deriv[:,1],Delta_eigvals_deriv[:,1]))
    
    fit_pi  = main_fitting(eig1,Nexp,T,tmin_pi,tmax_pi)
    fit_eta = main_fitting(eig2,Nexp,T,tmin_eta,tmax_eta)
    
    mpi  = fit_pi.p['dE'][0]
    meta = fit_eta.p['dE'][0]
    chi2perdof_pi  = fit_pi.chi2/fit_pi.dof
    chi2perdof_eta = fit_eta.chi2/fit_eta.dof 

    # Write everything into a csv 
    outHR = open(outfileHR, "a")
    outHR.write("%s;%s;%s;%s;%s;%s;%s;%s;%s;%s\n" % (h5file,Nexp,tmin_pi,tmax_pi,tmin_eta,tmax_eta,mpi,meta,chi2perdof_pi,chi2perdof_eta))
    outHR.close()

    out = open(outfile, "a")
    out.write("%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s\n" % (h5file,Nexp,tmin_pi,tmax_pi,tmin_eta,tmax_eta,gv.mean(mpi),gv.sdev(mpi),gv.mean(meta),gv.sdev(meta),chi2perdof_pi,chi2perdof_eta))
    out.close()


path = './output/eigvals/'
os.makedirs("output/data_human_readable/", exist_ok=True)
os.makedirs("output/data/", exist_ok=True)

outfileHR_name = "output/data_human_readable/non_degenerate_data_HR_eta_pi.csv"
outfile_name = "output/data/non_degenerate_data_eta_pi.csv"

outfileHR = open(outfileHR_name, "w")
outfileHR.write("file;Nexp;tmin;tmax;m_pi;m_eta';chi2perdof_pi;chi2perdof_eta\n")
outfileHR.close()

outfile = open(outfile_name, "w")
outfile.write("file;Nexp;tmin;tmax;m_pi;Delta_m_pi;m_eta';Delta_m_eta';chi2perdof_pi;chi2perdof_eta\n")
outfile.close()

with open('input/parameters/param_non_deg.csv') as csvfile:
    reader = csv.DictReader(csvfile,delimiter=';')
    for row in reader:
        # check if the file exists
        h5file  = path+os.path.basename(row['file'])+'.h5'
        if not os.path.isfile(h5file):
            continue        
        
        if row['fitπ0'] == "-1":
            continue

        tmin_pi,  tmax_pi  = literal_eval(row['fitπ0'])
        tmin_eta, tmax_eta = literal_eval(row['fitη']) 

        # perform fit
        Nexp = 5
        main(outfileHR_name,outfile_name,h5file,Nexp,tmin_pi,tmax_pi,tmin_eta,tmax_eta)