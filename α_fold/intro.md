1. [AlphaFold and Implications for Intrinsically Disordered Proteins](https://scholar.google.com/scholar?hl=en&as_sdt=0%2C5&q=AlphaFold+and+Implications+for+Intrinsically+Disordered+Proteins+Kiersten+M.+Ruff+and+Rohit+V.+pappus&btnG=&oq=AlphaFold+and+Implications+for+Intrinsically+Disordered+Proteins+Kiersten+M.+Ruff+and+Rohit+V.+Pappu)

Accurate predictions of the three-dimensional structures of proteins from their amino acid sequences have
come of age. AlphaFold, a deep learning-based approach to protein structure prediction, shows remarkable success in independent assessments of prediction accuracy. A significant epoch in structural bioinformatics was the structural annotation of over 98% of protein sequences in the human proteome.
Interestingly, many predictions feature regions of very low confidence, and these regions largely overlap
with intrinsically disordered regions (IDRs). That over 30% of regions within the proteome are disordered
is congruent with estimates that have been made over the past two decades, as intense efforts have been
undertaken to generalize the structure–function paradigm to include the importance of conformational
heterogeneity and dynamics. With structural annotations from AlphaFold in hand, there is the temptation
to draw inferences regarding the “structures” of IDRs and their interactomes.

Two recent reports highlight the coming of age of
machine learning and artificial intelligence in
structural bioinformatics.1,2
This success, as
assessed by adjudicators of the 14th edition of the
Critical Assessment of Protein Structure (CASP)
experiment, created quite a stir when it was
deployed to predict the structures of 98.5% of proteins in the human proteome. Through a portal
hosted by the European Bioinformatics Institute
(http://alphafold.ebi.ac.uk) one can now access predictions of protein structures made using AlphaFold. In an instant, detailed, atomic-level models
of protein structures have become available for
most of the protein sequences in the human proteome.

Interestingly, the AlphaFold predictions highlight
the importance of intrinsically disordered proteins/
regions (IDPs/IDRs). Conservative estimates
indicate that roughly 30% of sequences, 30-
residues or longer, drawn at random from the human proteome, are likely to be IDRs.3,5 A striking
feature of structural annotations of the human proteome provided by AlphaFold is the vast number
of low and very low confidence regions that overlap
with regions that are predicted to be IDRs. This
highlights the importance of conformational heterogeneity, which was advocated for in numerous
pioneering studies.6–8

Under standard conditions
used to study and characterize protein structures,
IDPs/IDRs are defined by conformational heterogeneity, and the preference for heterogeneity is
encoded in their amino acid sequences.11,29–31

Systematic efforts, spanning the past
two decades, have led to prediction engines that
allow one to identify regions within a protein
sequence that are likely to be IDRs.5,16,41–44
Given the advances enabled by AlphaFold, it is
now likely obligatory that biologists and
biochemists will look up the structures of their
favorite proteins. Many of these proteins,
especially
those
involved
in
signaling,
transcription, and coordinating protein–protein
interaction networks, are likely to feature large,
disordered regions. A typical annotation, shown in
Figure 1, will depict large regions as being orange
“unstructured” regions of very low confidence. The
key questions are: (1) What does one do with this
information? And (2) how should one interpret the
“unstructured regions” depicted by the AlphaFold
annotation?

One of the major
insights leveraged by AlphaFold1,2 comes from the
use of evolutionary covariations that can be
extracted from large-scale multiple sequence alignments.49–54 These discoveries, which go by various
names including direct coupling analysis (DCA),55
have their origins in early work showing that covariation analysis helps with the identification of functionally relevant sectors in protein structures.56
This is important to acknowledge, because one of
the challenges posed by IDRs stems from hyper-
variability of these regions across orthologs,57 making it difficult to uncover evolutionary constraints
from alignments alone.58

This is a per-residue confidence score
that is scaled between 0 and 100 and estimates
how well the predicted structure would agree with
the experimental structure as defined by the predicted Local Distance Difference Test (pLDDT).1,2
Here, we attempt to answer the following questions:
(1) Are all low and very low confidence regions failures of AlphaFold? i.e., do these regions adopt a
well-defined three-dimensional structure, but AlphaFold fails in predicting the relevant structures? (2)
Do the large clouds/ribbon-like depictions of very
low confidence regions in AlphaFold have physical
significance in terms of defining the radius of capture of the IDR or the interplay with folded domains?

Roughly 30% of the residues across predicted
structures in the human proteome tend to have
pLDDT scores that are less than 50.

it is
worth emphasizing that intrinsic disorder is not the
same as being “unstructured”. Instead, disorder
implies that a diverse conformational ensemble
best describes the region of interest. This
ensemble is sequence-specific30,31 and the extent
of heterogeneity36 as well as the relative preferences of conformations within the ensemble will
depend on the primary sequence, solution conditions, and functional contexts. The sequencespecificity of sequence-ensemble relationships
cannot be ignored. Accordingly, a single static
“structure”, even if it is annotated as being a low
confidence prediction, cannot be used as a representative conformation that describes the ensemble. Instead, what we need, and is being actively
pursued in the IDP field, are quantitative descriptions of conformational ensembles in terms of distribution functions for inter-residue distances and
measurements of moments of these distributions.61,62

Indeed, it is worth emphasizing that the
AlphaFold developers are explicit in making this
case, stating that: “In the current dataset, long
regions with pLDDT < 50 adopt a readily identifiable
ribbon-like appearance, and should not be interpreted as structures but rather as a prediction of
disorder.”

it is worth
noting that the conformational clouds one observes
in AlphaFold predictions may arise from several
aspects of the underlying methodology.1 First, multiple sequence alignments are needed to predict
distances between residues. Intrinsically disordered
regions often evolve more rapidly than ordered
regions and thus alignments of these regions are
generally poorer and involve large gaps and extensions of gaps because orthologous IDRs can span
vast sequence lengths.9,23,58,60,68–70 This likely
leads to poorly defined distance restraints within
the IDRs and between the IDRs and the folded
domains.

![fold1.png](https://github.com/LoqmanSamani/protein_sa/blob/systembiology/images/fold1.png)
![fols2.png](https://github.com/LoqmanSamani/protein_sa/blob/systembiology/images/fold2.png)

