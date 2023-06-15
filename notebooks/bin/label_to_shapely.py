#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2021 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""

"""
from shapely.geometry import Polygon, MultiPolygon
import cv2
import numpy as np
from scipy import ndimage
import tifffile as tf
import pickle
import fire


def get_shapely(label):
    """
    get outlines of masks as a list to loop over for plotting
    """
    polygons = {}
    simpler_polys = {}
    slices = ndimage.find_objects(label)
    for i, bbox in enumerate(slices):
        if not bbox:
            continue
        cur_cell_label = i + 1
        msk = (label[bbox[0], bbox[1]] == cur_cell_label).astype(np.uint8).copy()
        cnts, _ = cv2.findContours(
            msk, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE
        )
        # if len(cnts) > 1:
            # print(len(cnts), cur_cell_label)
        current_polygons = [
            Polygon((cnt + [bbox[1].start, bbox[0].start]).squeeze())
            #             else Point((cnt + [bbox[1].start, bbox[0].start]).squeeze())
            for cnt in cnts
            if len(cnt) > 2
        ]
        multipoly_obj = MultiPolygon(current_polygons)
        polygons[cur_cell_label] = multipoly_obj
        simpler_polys[cur_cell_label] = multipoly_obj.simplify(
            2, preserve_topology=True
        )
    return polygons, simpler_polys


def main(stem, label):
    lab = tf.imread(label).squeeze()
    print(lab.shape)
    with open(f"{stem}_cell_shapely.pickle", "wb") as handle:
        pickle.dump(
            get_shapely(lab)[0], handle, protocol=pickle.HIGHEST_PROTOCOL
        )


if __name__ == "__main__":
    fire.Fire(main)
