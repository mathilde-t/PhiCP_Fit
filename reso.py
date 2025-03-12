import awkward as ak

def compute_resolutions(events):
    hcands = [("hcand1", "h1", 1), ("hcand2", "h2", 2)]
    variables = ["pt", "eta", "phi"]
    comparisons = {
        "gen_reco": ("gen_{var}{num}", "{var}{num}"),
        "gen_fastmtt": ("gen_{var}{num}", "p4_{h_var}_reg.{var}"),
        "reco_fastmtt": ("p4_{h_var}_reg.{var}", "{var}{num}")
    }
    
    resolutions = {}
    for hcand, h_var, num in hcands:
        for var in variables:
            for comp, (num_fmt, denom_fmt) in comparisons.items():
                num_key = num_fmt.format(var=var, num=num, h_var=h_var)
                denom_key = denom_fmt.format(var=var, num=num, h_var=h_var)
                
                num_val = events[num_key] if num_key in events.fields else None
                denom_val = events[denom_key] if denom_key in events.fields else None
                
                if num_val is not None and denom_val is not None:
                    key = f"resolution_{hcand}_{var}_fastMTT_resolution_{comp}"
                    resolutions[key] = (num_val - denom_val) / denom_val if var == "pt" else num_val - denom_val
    
    # Concatenate results
    for comp in comparisons.keys():
        for var in variables:
            key_list = [resolutions[f"resolution_{hcand}_{var}_fastMTT_resolution_{comp}"] for hcand, _, _ in hcands]
            concatenated = ak.concatenate(key_list, axis=1)
            events = set_ak_column(events, f"hcand.{var}_fastMTT_W_resolution_{comp}", concatenated)
    
    return events
