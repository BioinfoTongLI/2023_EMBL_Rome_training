#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2020 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""
Assign peaks in .csv to labelled image
"""
import argparse
import pandas as pd
import pickle
import fire


def assign(trees, cells, stem="out"):

    spot_counts = {}
    cell_centroids = {}
    ys, xs, chs, cell_indexes = [], [], [], []
    for cell_index in cells:
        cell = cells[cell_index]
        if cell.centroid.is_empty: # ignore the shapes that don't have centoird
           continue
        cell_centroids[cell_index] = {"y": cell.centroid.y, "x": cell.centroid.x}

        current_counts = {}
        for ch in trees:
            potential_inside = trees[ch].query(cell)
            true_in = [trees[ch].geometries[i] for i in potential_inside if cell.is_valid and cell.contains(trees[ch].geometries[i])]
            for sp in true_in:
                ys.append(sp.y)
                xs.append(sp.x)
                chs.append(ch)
                cell_indexes.append(cell_index)
            current_counts[ch] = len(true_in)
        spot_counts[cell_index] = current_counts
        del cell
        del current_counts

    spots_df = pd.DataFrame(
        {"y": ys, "x": xs, "ch": chs, "ID": cell_indexes}
    ).set_index("ID")
    spots_df.to_csv(f"{stem}_assigned_peaks.csv")

    count_df = pd.DataFrame(spot_counts).T
    count_df.to_csv(f"{stem}_peak_counts.csv")

    centroid_df = pd.DataFrame(cell_centroids).T
    centroid_df.to_csv(f"{stem}_cell_centroids.csv")

    return count_df, centroid_df

    # n_total = [trees[i]._n_geoms for i in trees]
    # summary = pd.DataFrame(
    #     {
    #         "Sum": count_df.sum(),
    #         "Total_per_ch": n_total,
    #         "Percentage": count_df.sum() / n_total,
    #     }
    # )
    # summary.to_csv(f"{stem}_summary.csv")


if __name__ == "__main__":
    fire.Fire(assign)
