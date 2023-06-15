#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2021 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""
Get the STRtree of spots
"""
import argparse
import pandas as pd
from shapely.geometry import Point
from shapely.strtree import STRtree
import pickle
import fire


def get_STRtree_per_channel(spots_df, ch_col_name, x_col="x_int", y_col="y_col"):
    trees = {}
    spots_df.columns = spots_df.columns.str.lower()
    if ch_col_name != "":
        for i, spots in spots_df.groupby(ch_col_name):
            points = [
                Point(spots.loc[ind, x_col], spots.loc[ind, y_col])
                for ind in spots.index
            ]
            trees[i] = STRtree(points)
    else:
        points = [
            Point(spots_df.loc[ind, x_col], spots_df.loc[ind, y_col])
            for ind in spots_df.index
        ]
        trees["anchor"] = STRtree(points)
    return trees


def main(stem, peak, target_ch, sep, x_col, y_col):
    spots_df = pd.read_csv(peak, sep=str(sep))
    spots_df = spots_df[
        (spots_df[target_ch] != "background")
        & (spots_df[target_ch] != "infeasible")
        & (~spots_df[target_ch].isna())
        & (spots_df[target_ch] != "1.0")
    ]
    print(spots_df.shape)
    with open(f"{stem}_str_peaks.pickle", "wb") as handle:
        pickle.dump(
            get_STRtree_per_channel(spots_df, target_ch.lower(), x_col, y_col),
            handle,
            protocol=pickle.HIGHEST_PROTOCOL,
        )


if __name__ == "__main__":
    fire.Fire(main)
