# CP Hypothesis Fit and Asymmetry Analysis in Periodic Distributions

## Overview
A repository for plotting in the **CPinHToTauTau analysis**. This script performs a fit for three **CP hypotheses : <span style="color: #d62839;">**CP-even**</span>, <span style="color: #28348e;">**CP-odd**</span>, and <span style="color: #2b663c;">**maximal mixing**</span>. It the asymmetries between these hypotheses and visualizes the results are computed. The approach can be applied to a **general case where periodic distributions are fitted `OUTPUT/version_name/<process_name>/*` and their asymmetries are analysed `OUTPUT/Asymmetry_*.pdf`**.

## Features
- **Fit for CP Hypotheses :** The script provides fits to the data for the three CP hypotheses : CP-even, CP-odd, and maximal mixing.
- **Asymmetry Computation :** Asymmetries between these hypotheses for each output plot are computed.
- **Spin Analysing Power Visualization :** The results are plotted separately to analyse the spin-analyzing power of the underlying reconstruction method.

## Uses
Output (.pickle) from Columnflow

> ### How to prepare the files ?
> 1. **Copy the needed pickle files**
>    Use `scripts/phi_cp_merges_hists.sh` to:
>    - Select files by **CF production version**, **era**, and **process**.
>    - Copy them to your **CERNBox**.
>
> 2. **Organise the files locally**
>    Place the following structure inside the `INPUT` directory:
>
>    ```text
>    version_name/
>    ├── ggF
>    │   ├── sm
>    │   ├── mm
>    │   └── cpo
>    └── VBF
>        ├── sm
>        ├── mm
>        └── cpo
>    ```
>
> 3. **Rename the files**
>    Use `scripts/copy_merge_and_rename.sh` to:
>    - Remove **process** and **CP hypothesis** from filenames.
>    - => Ready to use for the **`PhiCPfit_BarPlot_DESYTau_cf_0p3.ipynb`** plotting script.


## Produces
- **Cosine fit** for Φ<sub>CP</sub> distributions in the **H → ττ → τ<sub>l</sub>τ<sub>h</sub>** channel for multiple hypothesis
- **Plots the Asymmetry** between different hypotheses

## Formulas for the Fitting and their Errors

1) ### Fit Formula  
   The **model function** for the fit analysis is :

   $`fit = f(x) = a \cdot \cos(x + c) + b`$

   Where $`a`$ is the amplitude, $`b`$ is the offset, and $`c`$ is the phase shift parameter.
  
2) ### Fit Parameter Errors  
   The uncertainties on the fit parameters $`a`$, $`b`$, and $`c`$ are determined using the **Hesse matrix** from the Minuit minimisation :

   $`\sigma_a, \sigma_b, \sigma_c = \text{Minuit.Hesse}(\text{fit})`$

   These errors are propagated into the asymmetry calculation.

3) ### Asymmetry Formula  
   The **asymmetry** between two distributions $`\text{fit}_{1}`$ and $`\text{fit}_{2}`$ is computed by :

   $`A_{1,2} = \frac{1}{N} \sum_{i=1}^{N} \left| \frac{\text{fit}_{1}^i - \text{fit}_{2}^i}{\text{fit}_{1}^i + \text{fit}_{2}^i} \right|`$

   Where $`N`$ is the number of elements in the arrays.

4) ### Asymmetry Error Formula  
   The error of the asymmetry is calculated using error propagation for the two distributions $`\text{fit}_{1}`$ and $`\text{fit}_{2}`$ :

   $`\sigma_A = \frac{1}{N} \sqrt{\sum_{i=1}^{N} \left( \left( \frac{2 \cdot \text{fit}_{2}^i \cdot \text{err}_1^i}{(\text{fit}_{1}^i + \text{fit}_{2}^i)^2} \right)^2 + \left( \frac{2 \cdot \text{fit}_{1}^i \cdot \text{err}_2^i}{(\text{fit}_{1}^i + \text{fit}_{2}^i)^2} \right)^2 \right)}`$

## Thanks 
Thanks to **Gourab Saha** for providing the **original Φ<sub>CP</sub> fit code** for the HToTauTau Analysis :  
[Gourab's Repo](https://github.com/gsaha009/PlayingCP)
