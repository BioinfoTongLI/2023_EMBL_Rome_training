#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2021 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""

"""
import pandas as pd
import re
import numpy as np
from collections import OrderedDict
import os


# nucleotide_map = {"1": "A", "2": "G", "3": "C", "4": "T"}


# def convert2AGCT(codelist):
    # AGCT_codes = []
    # for g in codelist:
        # AGCT_codes.append("".join([nucleotide_map[i] for i in str(g)]))
    # return AGCT_codes


# def get_channel_info(col_list):
    # n_chs = []
    # n_cycs = []
    # ch_info = {}
    # for col in col_list:
        # m = re.search("cycle(\d+)_channel(\d+)_(.*)", col)
        # if m:
            # cyc = int(m.group(1))
            # ch = int(m.group(2))
            # n_chs.append(ch)
            # n_cycs.append(cyc)
            # ch_name = m.group(3)
            # if ch_name == "DAPI":
                # ch_info[ch_name] = "nuclei"
            # else:
                # ch_info[ch_name] = nucleotide_map[str(ch)]
            # print(cyc, ch, ch_name, col)
    # n_ch_set = set(n_chs)
    # assert len(n_ch_set) == max(n_ch_set)
    # n_cyc_set = set(n_cycs)
    # assert len(n_cyc_set) == max(n_cyc_set)
    # ch_info["nCycles"] = max(n_cyc_set)
    # ch_info["nChannel"] = max(n_ch_set)
    # return ch_info


def main(csv_file,
         channel_map:dict={"Cy5": "A", "AF488": "G", "Cy3": "C", "Atto425": "T", "AF750":"T"}, # This is orderd! be cautious !!!
         sep=",",
         out_dir="./out/"):
    channel_map = OrderedDict(channel_map)
    print(channel_map)
    if csv_file.endswith(".xlsx"):
        d = pd.read_excel(csv_file)
    else:
        d = pd.read_csv(csv_file, sep=sep)
    code_sizes = [len(str(c)) for c in d.code]
    n_cycle_list = np.unique(code_sizes)
    assert len(n_cycle_list) == 1
    channel_info = {}
    channel_info["nCycles"] = n_cycle_list[0]

    channel_dict = {}
    for col in d.columns:
        if col.startswith("cycle"):
            m = re.search("cycle(\d+)_channel(\d+)_(.*)", col)
            channel_dict[(m.group(1), m.group(2))] = m.group(3)
    nucleotide_codes = []
    for gene_ind in d.index:
        gene = d.loc[gene_ind]
        str_l = []
        for i, ind in enumerate(str(gene.code)):
            ch_name = channel_dict[(str(i + 1), ind)]
            col_name = f"cycle{i+1}_channel{ind}_{ch_name}"
            assert gene[col_name] == 1
            str_l.append(ch_name)
        nucleotids = "".join([channel_map[s] for s in str_l])
        nucleotide_codes.append(nucleotids)
    d["nucleotide_codes"] = nucleotide_codes
    channel_indexes = [int(k[1]) for k in channel_dict.keys()]
    assert len(np.unique(channel_indexes)) == np.max(channel_indexes)
    channel_info["nChannel"] = np.max(channel_indexes)
    channel_info["DAPI"] = "nuclei"
    for ch in channel_map:
        if ch == "AF750":
            channel_info["Atto425"] = channel_map[ch]
        else:
            channel_info[ch] = channel_map[ch]
    # print(channel_info)
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    df = pd.DataFrame(pd.Series(channel_info)).T
    df.to_csv(f"{out_dir}/channel_info.csv", index=False)
    print(df)

    taglist = d[["gene", "nucleotide_codes"]]
    taglist = taglist.rename(columns={"gene": "Gene", "nucleotide_codes": "Channel"})
    print(taglist)
    pd.DataFrame(taglist).to_csv(f"{out_dir}/taglist.csv", index=False)

    # -------------------------------------------------------------------------
    # ch_info = get_channel_info(d.columns)
    # print(ch_info)
    # df = pd.DataFrame(pd.Series(ch_info)).T
    # columns_order = ["nCycles", "nChannel", "DAPI"]
    # nucleotides = [col for col in df.columns.values if col not in columns_order]
    # ordered_cols = columns_order + nucleotides
    # df = df[ordered_cols]

    # df.to_csv("channel_info.csv", index=False)
    # d["Channel"] = convert2AGCT(d.code)


if __name__ == "__main__":
    import fire
    fire.Fire(main)
