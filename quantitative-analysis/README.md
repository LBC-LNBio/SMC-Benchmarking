# Quantitative analysis

Here, we calculate the cavity volume with pyKVFinder for supramolecular cages B1 to B14.

## Commands

To execute pyKVFinder for the supramolecular cages, run:

```bash
pyKVFinder ../hosts/B1.pdb --step 0.25 --probe_out 10.0 --removal_distance 2.0 --volume_cutoff 5.0 --depth -O B1
pyKVFinder ../hosts/B2.pdb --step 0.25 --probe_out 10.0 --removal_distance 1.75 --volume_cutoff 5.0 --depth -O B2
pyKVFinder ../hosts/B3.pdb --step 0.25 --probe_out 10.0 --removal_distance 2.0 --volume_cutoff 5.0 --depth -O B3
pyKVFinder ../hosts/B4.pdb --step 0.25 --probe_out 10.0 --removal_distance 2.0 --volume_cutoff 5.0 --depth -O B4
pyKVFinder ../hosts/B5.pdb --step 0.25 --probe_out 10.0 --removal_distance 1.5 --volume_cutoff 5.0 --depth -O B5
pyKVFinder ../hosts/B6.pdb --step 0.25 --probe_out 10.0 --removal_distance 1.5 --volume_cutoff 5.0 --depth -O B6
pyKVFinder ../hosts/B7.pdb --step 0.25 --probe_out 10.0 --removal_distance 1.5 --volume_cutoff 110.0 --depth -O B7
pyKVFinder ../hosts/B8.pdb --step 0.25 --probe_out 10.0 --removal_distance 1.5 --volume_cutoff 25.0 --depth -O B8
pyKVFinder ../hosts/B9.pdb --step 0.25 --probe_out 10.0 --removal_distance 1.5 --volume_cutoff 5.0 --depth -O B9
pyKVFinder ../hosts/B10.pdb --step 0.25 --probe_out 10.0 --removal_distance 1.5 --volume_cutoff 5.0 --depth -O B10
pyKVFinder ../hosts/B11.pdb --step 0.25 --probe_out 10.0 --removal_distance 2.0 --volume_cutoff 5.0 --depth -O B11
pyKVFinder ../hosts/B12.pdb --step 0.25 --probe_out 10.0 --removal_distance 1.5 --volume_cutoff 5.0 --depth -O B12
pyKVFinder ../hosts/B13.pdb --step 0.25 --probe_out 6.0 --removal_distance 1.5 --volume_cutoff 5.0 --depth -O B13
pyKVFinder ../hosts/B14.pdb --step 0.25 --probe_out 10.0 --removal_distance 1.5 --volume_cutoff 5.0 --depth -O B14
```

## Cavities

To visualize cavities, load B.pse in PyMOL.

```bash
pymol B.pse
```

## Guest volumes


## Cavity volumes

To visualize cavity volumes, load B.xlsx.
