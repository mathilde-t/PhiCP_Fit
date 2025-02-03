import os
import re
#import ROOT
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
import itertools

tags = ["IPIP", "IPDP", "IPPV"]  # Define the tags
catdict = {
    r"$\mu-\pi$"     : 20014002,
    r"$\mu-\rho$"    : 20014004,
    r"$\mu-a^{1}_{3\pi}$" : 20014008,
}
shiftdict = {
    "cp_even": {"shift": 150, "colour": "red", "location": "upper right"},
    "cp_odd": {"shift": 151, "colour": "blue", "location": "lower left"},
    "cp_maxmix": {"shift": 0, "colour": "black", "location": "lower right"},
}
simpledict = {
    r"$\mu-\pi$"     : "mupi",
    r"$\mu-\rho$"    : "murho",
    r"$\mu-a^{1}_{3\pi}$" : "mua13pr",
}

fit_results = { "tags": {}, "cats": {}, "asymmetry": {} }         # e.g tag = IPDP, cat == key = mupi, Asymmetry = num
asymmetry_results = {}  # Dictionary to store asymmetry results per tag

def comp_asymmetry(arr1, arr2):
    # https://github.com/Ksavva1021/TIDAL/blob/656f992ae056b3fed0061f2b3efb49905c39834d/CP_Tools/scripts/assymetry.py#L26
    return (1/arr1.size)*np.sum(np.abs((arr1-arr2)/(arr1+arr2)))

def comp_asymmetry_error(arr1, arr2, err1, err2): #arr = array of values, err = array of errors
    # Terms for partial derivatives
    denom = arr1 + arr2
    term1 = err1 * np.abs((2 * arr2) / (denom**2))
    term2 = err2 * np.abs((2 * arr1) / (denom**2))
    # Propagate errors
    sigma_A = np.sqrt(np.sum(term1**2 + term2**2)) / arr1.size
    return sigma_A

for tag in tags:
    # Update the file path according to the tag
    file = f"INPUT/shifted_hist__PhiCPGen_{tag}.pickle"

    # Load data and perform calculations as in your original code (already done in your script)
    fileptr = open(file, 'rb')
    data = pickle.load(fileptr)
    fileptr.close()

    axes = data.axes
    category_axis = axes['category']
    shift_axis = axes['shift']
    
    cparray = {}

    for ccat, cval in catdict.items():
        shiftarray = {}
        for cat, props in shiftdict.items():
            shift = props["shift"]
            colour = props["colour"]
            location = props["location"]

            if cval not in category_axis:
                print(f"WARNING : {cval} not in categories")
                continue

            shift_index = shift_axis.index(props["shift"])    
            values = data[category_axis.index(cval), :, shift_axis.index(props["shift"]), :].values()
            errors = data[category_axis.index(cval), :, shift_axis.index(props["shift"]), :].variances() ** 0.5
            shiftarray[cat] = {
                "values": values,
                "errors": errors,
                "colour": colour,
                "location": location,
            }

        cparray[ccat] = shiftarray

    # Now calculate the asymmetry values and store them in the results dictionary for the current tag
    for category1, category2 in itertools.combinations(shiftdict.keys(), 2): # Loop over pairs of categories
        hypothesis1 = np.ravel(cparray[ccat][category1]["values"])
        hypothesis2 = np.ravel(cparray[ccat][category2]["values"])
        error1 = np.ravel(cparray[ccat][category1]["errors"])
        error2 = np.ravel(cparray[ccat][category2]["errors"])

        asymmetry = comp_asymmetry(hypothesis1, hypothesis2)
        asymmetry_error = comp_asymmetry_error(hypothesis1, hypothesis2, error1, error2)

        # Store the results in the dictionary under the tag and category pair
        asymmetry_results[tag] = asymmetry_results.get(tag, {})
        asymmetry_results[tag][f"{category1}_vs_{category2}"] = {
            "asymmetry_val": asymmetry,
            "asymmetry_error": asymmetry_error
        }

for tag in tags:
    plt.figure(figsize=(8.9, 6.6))
    hep.cms.text("Simulation", loc=1)

    # Prepare the data for plotting
    asymmetry_vals = []
    asymmetry_errs = []
    categories = []

    for category in catdict.keys():  # Loop through all categories
        # Retrieve the asymmetry value and error for each category under the current tag
        asymmetry_val = asymmetry_results[tag].get(f"cp_even_vs_cp_odd", {}).get("asymmetry_val", 0)
        asymmetry_error = asymmetry_results[tag].get(f"cp_even_vs_cp_odd", {}).get("asymmetry_error", 0)

        if asymmetry_val != 0:
            asymmetry_vals.append(asymmetry_val)
            asymmetry_errs.append(asymmetry_error)
            categories.append(simpledict.get(category, category))  # Use 'simpledict' to convert category to label

    # Plot the results for this tag
    plt.errorbar(categories, asymmetry_vals, yerr=asymmetry_errs, fmt="o", label=f"Tag: {tag}")

    # Customize the plot
    plt.xlabel("Categories", fontsize=20)
    plt.ylabel("Asymmetry Value", fontsize=20)
    plt.title(f"Asymmetry for {tag}", fontsize=25)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.legend()

    # Save the plot
    plt.savefig(f"OUTPUT/{tag}_asymmetry_plot.pdf", dpi=300)
    plt.show()
