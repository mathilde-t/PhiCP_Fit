import os
import re
import ROOT
import hist
import pickle
import mplhep as hep
from glob import glob
import numpy as np
import awkward as ak
from tabulate import tabulate
from IPython import embed
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
# iminuit
import iminuit
from iminuit import Minuit
from iminuit.cost import LeastSquares
hep.style.use("CMS")
print("iminuit version:", iminuit.__version__)

tag = "IPNP" #q? what is this for ?

catdict = {
    r"$\mu-\pi$"     : 401311, #q? can i leave the tag numbers in here ?
    r"$\mu-\rho$"    : 401313,
    r"$\mu-a^{1}_{3\pi}$" : 401317,

}
shiftdict = {
    "cp_even" : 150, #q? can i leave the tag numbers in here ?
    "cp_odd"  : 151,
}
simpledict = {
    r"$\mu-\pi$"     : "mupi",
    r"$\mu-\rho$"    : "murho",
    r"$\mu-a^{1}_{3\pi}$" : "mua13pr",
}

def comp_asymmetry(arr1, arr2):
    # https://github.com/Ksavva1021/TIDAL/blob/656f992ae056b3fed0061f2b3efb49905c39834d/CP_Tools/scripts/assymetry.py#L26
    return (1/arr1.size)*np.sum(np.abs((arr1-arr2)/(arr1+arr2)))

def makesimple(latex_str):
    # Remove LaTeX commands (e.g., $...$, \frac, \text, etc.)
    plain_text = re.sub(r'\\[a-zA-Z]+\{[^}]*\}', '', latex_str)  # Remove \command{...}
    print(plain_text)
    plain_text = re.sub(r'\$.*?\$', '', plain_text)  # Remove math expressions between $...$
    print(plain_text)
    plain_text = re.sub(r'\\[a-zA-Z]+', '', plain_text)  # Remove other LaTeX commands (e.g., \chi)
    print(plain_text)

    # Clean up any remaining special characters
    plain_text = plain_text.replace('$', '')  # Remove leftover dollar signs
    print(plain_text)    
    return plain_text

np.linspace(0., 10., 20)

file = f"INPUT/shifted_hist__PhiCPGen_{tag}.pickle" #TODO change this file
#q? which file to use ? -> /eos/user/m/mawitt/share/OutputCP/cf_store/analysis_httcp/cf.PlotShiftedVariables1D/run3_2022_preEE_nano_cp_tau_v14/calib__main/sel__main/prod__main/weight__main/shifts_tauspinner/datasets_h_ggf_tautau_uncorrelated_filter/PhiCPtt+mt_12jan25


fileptr = open(file, 'rb')
data = pickle.load(fileptr)
fileptr.close()

axes = data.axes
category_axis  = axes['category']
shift_axis = axes['shift']

cparray = {}
for ckey, cval in catdict.items():
    shiftarray = {}
    for key,val in shiftdict.items():
        if cval not in category_axis:
            print(f"WARNING : {cval} not in categories")
            continue
        values = data[category_axis.index(cval), :, shift_axis.index(val), :].values()
        # https://github.com/oponcet/CPinHToTauTau/blob/FF_dev_project/script_FF/fake_factor_derivation/src/input_processing.py#L133
        errors = data[category_axis.index(cval), :, shift_axis.index(val), :].variances() ** 0.5
        #shiftarray[key] = data[category_axis.index(cval), :, shift_axis.index(val), :].values()
        shiftarray[key] = [values, errors]
    cparray[ckey] = shiftarray

cpfitarray = {}
x = np.linspace(0., 2*np.pi, 20)
def model(x, a, b):
    return a*np.cos(x) + b

def fit(x, y, err=0.05, model=model):
    lsq = LeastSquares(x, y, err, model)
    m = Minuit(lsq, a=0.1, b=0.1)
    #m.scan(ncall=100)
    m.fixed = False
    m.migrad()  # finds minimum of least_squares function
    m.hesse()  # accurately computes uncertainties    
    return m, err

for key, val in cparray.items():
    print(key)
    if len(val) == 0:
        print(f"WARNING : {key} has empty dict")
        continue
    even_zip = val["cp_even"]
    even, even_err = even_zip[0][0], even_zip[1][0]
    odd_zip  = val["cp_odd"]
    odd, odd_err = odd_zip[0][0], odd_zip[1][0]
    m_even, err_even = fit(x, even, even_err)
    m_odd, err_odd = fit(x, odd, odd_err)

    cpfitarray[key] = {"m_even": m_even, "m_odd": m_odd}

    # plot
    plt.figure(figsize=(8.9, 6.6))
    #plt.subplots_adjust(top=0.85)

    # Add CMS-style text
    hep.cms.text("Simulation", loc=1)
    #plt.text(0, 21, "CMS Simulation", fontsize=12, ha='left')
    
    plt.errorbar(x, even, even_err, fmt="o",color="blue")
    #plt.scatter(x, even, color="blue")
    even_fit = model(x, *m_even.values)
    lin1, = plt.plot(x, even_fit, color="blue")

    # display legend with some fit info
    fit_info_even = [
        f"$\\chi^2$/$n_\\mathrm{{dof}}$ = {m_even.fval:.1f} / {m_even.ndof:.0f} = {m_even.fmin.reduced_chi2:.1f}",
    ]

    #for p, v, e in zip(m.parameters, m.values, m.errors):
    #    fit_info.append(f"{p} = ${v:.3f} \\pm {e:.3f}$")
    leg_even_handle = Line2D([0], [0], color='blue', label="CP even")
    leg_even = plt.legend(handles=[leg_even_handle],title="\n".join(fit_info_even), frameon=False, loc="upper right", fontsize=20, title_fontsize=15)
    plt.gca().add_artist(leg_even)
    

    plt.errorbar(x, odd, odd_err, fmt="o",color="red")
    #plt.scatter(x, odd, color="red")
    odd_fit = model(x, *m_odd.values)
    lin2, = plt.plot(x, odd_fit, color="red")

    fit_info_odd = [
        f"$\\chi^2$/$n_\\mathrm{{dof}}$ = {m_odd.fval:.1f} / {m_odd.ndof:.0f} = {m_odd.fmin.reduced_chi2:.1f}",
    ]

    leg_odd_handle = Line2D([0], [0], color='red', label="CP odd")
    plt.legend(handles=[leg_odd_handle], title="\n".join(fit_info_odd), frameon=False, loc="lower left", fontsize=20, title_fontsize=15)
    #plt.gca().add_artist(leg_odd)

    plt.xlabel(r"$\Phi_{CP}$"+f" ({tag})")
    plt.ylabel("a.u")

    asymm = comp_asymmetry(even_fit, odd_fit)
    
    plt.title(f"{key} (A = {round(asymm, 3)})", fontsize=25, loc='center')
    

    plt.ylim(0.0, 
             np.max(odd_fit) + 0.65*np.max(odd_fit))

    plt.tight_layout()
    plt.savefig(f"OUTPUT/{tag}_{simpledict[key]}.pdf", dpi=300)
    plt.show()