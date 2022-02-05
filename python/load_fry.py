import scanpy

def load_fry(frydir, which_counts=['S','A'], verbose=False):
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
    elif verbose:
        print("processing input in standard mode, will return spliced count")

    af_raw = scanpy.read_mtx(os.path.sep.join([frydir, "alevin", "quants_mat.mtx"]))
    afg = [ l.rstrip() for l in open(os.path.sep.join([frydir, "alevin", "quants_mat_cols.txt"])).readlines()][:ng]
    afg_df =  pd.DataFrame(afg, columns=["gene_ids"])
    afg_df = afg_df.set_index("gene_ids")
    
    abc = [ l.rstrip() for l in open(os.path.sep.join([frydir, "alevin", "quants_mat_rows.txt"])).readlines() ]
    abc_df = pd.DataFrame(abc, columns=["barcodes"])
    abc_df.index = abc_df["barcodes"]
    
    x = af_raw.X
    if usa_mode:
        rd = {'S' : range(0,ng), 'U' : range(ng, 2*ng), 'A' : range(2*ng,3*ng)}
        o = x[:, rd[which_counts[0]]]
        for wc in which_counts[1:]:
            o += x[:, rd[wc]]
    else:
        o = x
        
    af = scanpy.AnnData(o.T, var=abc_df, obs=afg_df)
    af = af.T
    return af
