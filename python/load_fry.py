import scanpy

def load_fry(frydir, which_counts={'X' : ['S','A']}, verbose=False):
    """
    
    Parameters:
        frydir - The directory containing the alevin-fry quantification (i.e. the the quant.json file & alevin subdirectory).
        verbose - True if messages (including error messages) should be printed out, False if function should be quiet.
        which_count - Dictionary specifying how a USA mode matrix should be returned or combined into the resulting 
                      output matrix.  If the input is not a USA mode quantification directory, this parameter is ignored
                      and the count matrix is returned in the `X` field of the returned `AnnData` object.  If the input
                      quantification directory contains a USA mode quantification, then there are 3 sub-matrices that can 
                      be referenced in the dictionary; 'U', 'S', 'A' containing, respectively, unspliced, spliced and 
                      ambiguous counts.  The dictionary should have entries of the form `key` (str) : `value` (list[str]).
                      The following constraints apply : there should be one key-value pair with the key `X`, the resulting
                      counts will be returned in the `X` field of the AnnData object. There can be an arbitrary number
                      of other key-value pairs, but each will be returned as a layer of the resulting AnnData object.
                      Within the key-value pairs, the key refers to the layer name that will be given to the combined 
                      count matrix upon output, and the value should be a subset of `['U', 'S', 'A']` that defines 
                      which sub-matrices should be summed.  For example:
                      {'X' : ['S', 'A'], 'unspliced' : ['U']}

                      will result in a return AnnData object where the X field has a matrix in which each entry 
                      corresponds to the summed spliced and ambiguous counts for each gene in each cell, and there
                      is an additional 'unspliced' layer, whose counts are taken directly from the unspliced sub-matrix.

    Returns:
        An AnnData object with X and layers corresponding to the requested `which_counts`, or None if an 
        error is encountered.
    """
    import json
    import os
    import pandas as pd

    # since alevin-fry 0.4.1 the generic "meta_info.json"
    # has been replaced by a more informative name for each
    # sub-command. For quantification, it is "quant.json".
    # we check for both files here, in order.
    meta_info_files = ["quant.json", "meta_info.json"]

    fpath = os.path.sep.join([frydir, meta_info_files[0]])
    # first, check for the new file, if we don't find it, check
    # for the old one.
    if not os.path.exists(fpath):
        if verbose:
            print(f"Did not find a {meta_info_files[0]} file, checking for older {meta_info_files[1]}.")
        fpath = os.path.sep.join([frydir, meta_info_files[1]])
        # if we don't find the old one either, then return None
        if not os.path.exists(fpath):
            if verbose:
                print(f"Found no {meta_info_files[1]} file either; cannot proceed.")
            return None

    # if we got here then we had a valid json file, so 
    # use it to get the number of genes, and if we are 
    # in USA mode or not.
    meta_info = json.load(open(fpath))
    ng = meta_info['num_genes']
    usa_mode = meta_info['usa_mode']

    # if we are in USA mode
    if usa_mode:
        # make sure that num_genes is a multiple of 3
        if ng %3 != 0:
            if verbose:
                print("Found USA mode, but num genes = {ng} is not a multiple of 3; cannot proceed.")
            return None
        # each gene has 3 splicing statuses, so the actual number of distinct 
        # genes is ng/3.
        ng = int(ng/3)
        if verbose:
            print("processing input in USA mode, will return {}".format("+".join(which_counts)))
              
        # make sure which_counts isn't empty
        assert(len(which_counts) > 0)  

        # make sure the specification in which_counts is OK
        if 'X' not in which_counts:
            if verbose:
                print('In USA mode some sub-matrices must be assigned to the \"X\" (default) output.')
            return None
        if verbose:
            print(f"will populate output field X with sum of counts frorm {which_counts['X']}.")

        for k,v in which_counts.items():
            valid_elem = len(set(v) - set(['U', 'S', 'A'])) == 0
            if not valid_elem:
                if verbose:
                    print(f'Found non-USA element in which_count element list \"{v}\" for key \"{k}\"; cannot proceed.')
                return None
            if verbose:
                print(f'will combine {v} into output layer {k}.') 

    elif verbose:
        print("Processing input in standard mode, will return processed count (which_count will be ignored).")

    # read the actual input matrix
    af_raw = scanpy.read_mtx(os.path.sep.join([frydir, "alevin", "quants_mat.mtx"]))
    afg = [ l.rstrip() for l in open(os.path.sep.join([frydir, "alevin", "quants_mat_cols.txt"])).readlines()][:ng]
    # read the gene ids
    afg_df =  pd.DataFrame(afg, columns=["gene_ids"])
    afg_df = afg_df.set_index("gene_ids")
    # and the barcodes
    abc = [ l.rstrip() for l in open(os.path.sep.join([frydir, "alevin", "quants_mat_rows.txt"])).readlines() ]
    abc_df = pd.DataFrame(abc, columns=["barcodes"])
    abc_df.index = abc_df["barcodes"]
    
    x = af_raw.X
    # if we're not in USA mode, just combine this info into 
    # an AnnData object
    if not usa_mode:
        af = scanpy.AnnData(x.T, var=abc_df, obs=afg_df)
        af = af.T 
    else: # USA mode
        # otherwise, combine the sub-matrices into the output object as 
        # specified by `which_counts`
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
