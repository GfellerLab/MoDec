# Release notes from MoDec

## Version 1.2 (23.11.2022)

* Added optional arguments `MHC2`, `no_bestPepResp`, `kmer`, `alphabet`,
  `unk_aa`, `col_scheme`, `theta_norm`, `flat_freq`, `Salign` (see the README
  for details about these options).
* The number of peptides below each motif in the report now corresponds to
  number of peptides assigned to this motif, after assigning each peptide to
  one motif based on its highest responsibility value (instead of considering
  fractional values based on the peptides's weight and responsibilities).
* The flat motif now has the same length as the other motifs (it does not
  show anymore the values used for this flat motif in the deconvolution as
  it was confusing to see motifs of different sizes).
* Updated license-related contact person.

## Version 1.1 (14.10.2019)

* Improved the layout of the html report, of the readme and other small changes.

## Version 1.0 (11.02.2019)

* Public release of MoDec version 1.0.
