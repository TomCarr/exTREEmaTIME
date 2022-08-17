# exTREEmaTIME
A method for incorporating uncertainty into divergence time estimates.\
Instructions and further details are [**here**](https://github.com/TomCarr/exTREEmaTIME/wiki/Further-details-and-instructions).\
please cite: Carruthers and Scotland 2022. exTREEmaTIME: a method for incorporating uncertainty into divergence time estimates. Biol Open. https://doi.org/10.1242/bio.059181 

### 17/08/22
examples clarified. The original examples (v0) don't work with v1.0 (because it searches for the discontinued argument of calibration_implementation_precision). Example provided that should work straight away with version 1.

### update 20/12/21
Version 1 added. Runs about 1000X faster than the original version. Produces virtually identical estimates. calibration_implementation_precision parameter removed. All calibrations are implemented precisely.  

### update 24/11/21
Examples updated in accordance with below.

### update 15/11/21
We have become aware of one minor error in the constraints implemented in the empirical examples provided. Specifically, Amborellaceae has been missed from the angiosperm crown node maximum age constraint. This error has essentially no effect on divergence time estimates apart from the timing of the divergence between Amborella and the rest of Angiosperms. However, for clarity, updated analyses with the correct constraints will be uploaded later today.
