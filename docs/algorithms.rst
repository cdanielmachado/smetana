==========
Algorithms
==========

SMETANA implements several algorithms to analyse cross-feeding interactions in microbial communities. These have been
describe in `Zelezniak et al, PNAS (2015) <http://www.pnas.org/content/112/20/6449.short>`_. Please read the paper for
a more detailed explanation.


Global
______

SMETANA includes two algorithms that analyse global properties of a community:

- MRO (metabolic resource overlap): calculates how much the species compete for the same metabolites.
- MIP (metabolic interaction potential): calculates how many metabolites the species can share to decrease their dependency on external resources.
- SMETANA: the global SMETANA score is a measure of the community interactions and is the sum of the individual SMETANA scores.



Detailed
________

SMETANA includes several scores to characterize individual interactions in a community. This includes:

- SCS (species coupling score): measures the dependency of one species in the presence of the others to survive
- MUS (metabolite uptake score): measures how frequently a species needs to uptake a metabolite to survive
- MPS (metabolite production score): measures the ability of a species to produce a metabolite
- SMETANA: the individual smetana score is a combination of the 3 scores above, it gives a measure of certainty on a cross-feeding interaction (species A receives metabolite X from species B).

