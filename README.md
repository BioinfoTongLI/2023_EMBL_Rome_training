# 2023_EMBL_Rome_training
matreial for EMBO Rome image analysis training

__Disclaimer:__

This repo is only for educational purpose and will have troubles in production. Mainly because of the lack of scalability part.
And for each of the step, you will find more than one tool of doing it.
For example,
- [Stardist](https://qupath.readthedocs.io/en/0.3/docs/advanced/stardist.html) for cell segmentation
- [StarFish](https://spacetx-starfish.readthedocs.io/en/latest/) for tile-based preprocessing/peak-calling/decoding
- [TissUUMaps](https://tissuumaps.github.io/) for spatial transcriptomics data visualisation
- [Ashlar](https://github.com/labsyspharm/ashlar) for image stitching/registration
  
We select the ones in this repo because they are shown to be at the sweet point between precision and robustness with different tissues.


## Outline

1. Cycle-wise image regitration with [Microaligner](https://github.com/VasylVaskivskyi/microaligner)
2. Peak-calling with [TrackPy](http://soft-matter.github.io/trackpy/v0.6.1/)&Peak-decoding with [PoSTcode](https://github.com/milana-gataric/postcode)
3. Cell segmentation with [CellPose 2.0](https://www.cellpose.org/)
4. Assignment of detected barcodes to segmented cells with [STRtree](https://shapely.readthedocs.io/en/stable/strtree.html)
5. Data visualization with [Webatlas pipeline](https://developmental.cellatlas.io/webatlas)
