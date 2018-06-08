# Useful scripts for molecular biology applications

## primer_design.py

Automated primer design using Primer3 (can be used as a standalone script or loaded as a module)

## ugene_fetch.py

Creates a Ugene (also compatible with ApE) of a specified *C. elegans* genomic region, including gene features. Alternatively, add features from a gff/gtf file to a pre-existing Ugene file created using ugene_fetch.

### Examples:

Make Ugene file for trt-1
```
ugene_fetch.py -c I -s 8768325 -e 8771248 -n trt1
```
Add features to trt-1 file
```
ugene_fetch.py -g features.gtf -a trt1.gb
```
