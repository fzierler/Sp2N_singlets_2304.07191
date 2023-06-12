# Singlets in gauge theories with fundamental matter
This repository contains the code used to prepare the plots and results included in [Singlets in gauge theories with fundamental matter [2304.07191]](https://arxiv.org/abs/2304.07191v1).

## Instructions: Running the analysis
- Install required dependencies (see below)
- Download the hdf5-file `data.hdf5` from Zenodo and place it in `input_data` (the directory can be modified in `main.jl`)
- Download the file `WilsonFlow.zip`, decompress and place it in `input_data` (the directory can be modified in `main.jl`)
- Run the analysis using `julia main.jl` within this directory
- The figures and tables of [[2304.07191]](https://arxiv.org/abs/2304.07191v1) can then be found in 
    - `output/figures`
    - `output/tables`

- If you want to start from the raw logs the variable `start_from_raw_logs` in the file `main.jl` needs to be set to `true` and the path to the directory containing the decompressed raw logs needs to be provided. Note, that the raw logs are compressed. The variable `start_from_raw_logs` is set to `false` by default.

## Warning

The code in this repository has only been tested on the specific dataset provided here. It is not intended to be easily generalizable to arbitrary datasets. The analysis parameters are hard-coded in the directory `input/parameters`.

## Plots

The plots are made using [Plots.jl](https://zenodo.org/record/7994271) via the [PGFPlotsX](https://github.com/KristofferC/PGFPlotsX.jl) backend which requires a LaTeX installation with the PGFPlots package.

## Gradient flow analysis

The gradient flow data is analysed using Ed Bennett's [flow_analysis](https://github.com/edbennett/flow_analysis). 

## Requirements
- Python 3.8 (numpy, pre-commit, scipy, uncertainties)
- julia 1.9
- LaTeX (including PGFPlots)

## References for SU(3) 

In [[2304.07191]](https://arxiv.org/abs/2304.07191v1) we collected results from SU(3) gauge theory with two fundamental fermions. The data is given in tabulated form in `input/su3_literature.csv`.

In some cases, the measurement has been performed using different methods in the analysis or different operators have been used to study the same mesons and sets of results are available. In such cases, we have chosen the results that are closest to the determination of directly fitting the correlator of a pure fermionic operator. When this was not possible we quote the largest and smallest values of all measurements i and symmetrize the uncertainties. We detail our choices in `input/README_su3_literature.csv`

We stress that this is only meant to provide a qualitative understanding of the fermion-mass dependence of the pseudoscalar singlet meson mass. Our way of extracting the data might have introduced artificially large uncertainties and will certainly have introduced biases. Systematic uncertainties were only taken into account if the original reference provided them. For quantitative results and discussions of the systematic uncertainties we refer to the original literature.
