import scanpy

def load_fry(frydir, which_counts={'X' : ['S','A']}, verbose=False):
    import json
    import os
    import pandas as pd

    # since alevin-fry 0.4.1 the generic "meta_info.json"
    # has been replaced by a more informative name for each
    # sub-command. For quantification, it is "quant.json".
    # we check for both files here, in order.
    meta_info_files = ["quant.json", "meta_info.json"]

    fpath = os.path.sep.join([frydir, meta_info_files[0]])
    if not os.path.exists(fpath):
        if verbose:
            print(f"Did not find a {meta_info_files[0]} file, checking for older {meta_info_files[1]}.")
        fpath = os.path.sep.join([frydir, meta_info_files[1]])
        if not os.path.exists(fpath):
            if verbose:
                print(f"Found no {meta_info_files[1]} file either; cannot proceed.")
            return None

    meta_info = json.load(open(fpath))
    ng = meta_info['num_genes']
    usa_mode = meta_info['usa_mode']

    if usa_mode:
        assert(len(which_counts) > 0)
        if ng %3 != 0:
            if verbose:
                print("Found USA mode, but num genes = {ng} is not a multiple of 3; cannot proceed.")
            return None
        ng = int(ng/3)
        if verbose:
            print("processing input in USA mode, will return {}".format("+".join(which_counts)))
        
        # make sure the specification in which_counts is OK
        if 'X' not in which_counts:
            if verbose:
                print('In USA mode some sub-matrices must be assigned to the \"X\" (default) output.')
            return None
        for k,v in which_counts.items():
            valid_elem = len(set(v) - set(['U', 'S', 'A'])) == 0
            if not valid_elem:
                if verbose:
                    print(f'Found non-USA element in which_count element list \"{v}\" for key \"{k}\"; cannot proceed.')
                return None 

    elif verbose:
        print("Processing input in standard mode, will return processed count (which_count will be ignored).")

    af_raw = scanpy.read_mtx(os.path.sep.join([frydir, "alevin", "quants_mat.mtx"]))
    afg = [ l.rstrip() for l in open(os.path.sep.join([frydir, "alevin", "quants_mat_cols.txt"])).readlines()][:ng]
    afg_df =  pd.DataFrame(afg, columns=["gene_ids"])
    afg_df = afg_df.set_index("gene_ids")
    
    abc = [ l.rstrip() for l in open(os.path.sep.join([frydir, "alevin", "quants_mat_rows.txt"])).readlines() ]
    abc_df = pd.DataFrame(abc, columns=["barcodes"])
    abc_df.index = abc_df["barcodes"]
    
    x = af_raw.X
    if not usa_mode:
        af = scanpy.AnnData(x.T, var=abc_df, obs=afg_df)
        af = af.T 
    else: # USA mode
        rd = {'S' : range(0,ng), 'U' : range(ng, 2*ng), 'A' : range(2*ng,3*ng)}
        xcounts = which_counts['X']
        o = x[:, rd[xcounts[0]]]
        for wc in xcounts[1:]:
            o += x[:, rd[wc]]
        af = scanpy.AnnData(o.T, var=abc_df, obs=afg_df)
        af = af.T

        # now, if there are other layers requested, populate those
        for other_layer in which_counts.keys() - 'X':
            xcounts = which_counts[other_layer]
            o = x[:, rd[xcounts[0]]]
            for wc in xcounts[1:]:
                o += x[:, rd[wc]] 
            af.layers[other_layer] = o
    return af
