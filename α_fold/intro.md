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


2. [Highly accurate protein structure prediction for the human proteome](https://scholar.google.com/scholar?hl=en&as_sdt=0%2C5&q=Highly+accurate+protein+structure+prediction+for+the+human+proteome+Kathryn+Tunyasuvunakool1+%E2%9C%89%2C+Jonas+Adler1%2C+Zachary+Wu1%2C+Tim+Green1%2C+Michal+Zielinski1%2C+Augustin+%C5%BD%C3%ADdek1%2C+Alex+Bridgland1%2C+Andrew+Cowie1%2C+Clemens+Meyer1%2C+Agata+Laydon1%2C+Sameer+Velankar2%2C+Gerard+J.+Kleywegt2%2C+Alex+Bateman2%2C+Richard+Evans1%2C+Alexander+Pritzel1%2C+Michael+Figurnov1%2C+Olaf+Ronneberger1%2C+Russ+Bates1%2C+Simon+A.+A.+Kohl1%2C+Anna+Potapenko1%2C+Andrew+J.+Ballard1%2C+Bernardino+Romera-Paredes1%2C+Stanislav+Nikolov1%2C+Rishub+Jain1%2C+Ellen+Clancy1%2C+David+Reiman1%2C+Stig+Petersen1%2C+Andrew+W.+Senior1%2C+Koray+Kavukcuoglu1%2C+Ewan+Birney2%2C+Pushmeet+Kohli1%2C+John+Jumper1%2C3+%E2%9C%89+%26+Demis+Hassabis1%2C3+%E2%9C%89&btnG=)

Protein structures can provide invaluable information, both for reasoning about
biological processes and for enabling interventions such as structure-based drug
development or targeted mutagenesis. After decades of effort, 17% of the total
residues in human protein sequences are covered by an experimentally determined
structure1. Here we markedly expand the structural coverage of the proteome by
applying the state-of-the-art machine learning method, AlphaFold2, at a scale that
covers almost the entire human proteome (98.5% of human proteins). The resulting
dataset covers 58% of residues with a confident prediction, of which a subset (36% of
all residues) have very high confidence.

Thanks to the efforts of individual
laboratories and dedicated structural genomics initiatives, more than
50,000 human protein structures have now been deposited, making
Homo sapiens by far the best represented species in the Protein Data
Bank (PDB)5. Despite this intensive study, only 35% of human proteins
map to a PDB entry, and in many cases the structure covers only a
fragment of the sequence6. Experimental structure determination
requires overcoming many time-consuming hurdles: the protein must
be produced in sufficient quantities and purified, appropriate sample
preparation conditions chosen and high-quality datasets collected. A
target may prove intractable at any stage, and depending on the chosen
method, properties such as protein size, the presence of transmem-
brane regions, presence of disorder or susceptibility to conformational
change can be a hindrance7,8.

Protein structure prediction contributes to closing this gap by pro-
viding actionable structural hypotheses quickly and at scale. Previ-
ous large-scale structure prediction studies have addressed protein
families9–12, specific functional classes13,14, domains identified within
whole proteomes15 and, in some cases, full chains or complexes16,17.

particular, projects such as the SWISS-MODEL Repository, Genome3D
and ModBase have made valuable contributions by providing access
to large numbers of structures and encouraging their free use by the
community17–19. Related protein bioinformatics fields have developed
alongside structure prediction, including protein design20,21, function
annotation22–24, disorder prediction25, and domain identification and
classification26–28.

0. In particular, the latest version of
AlphaFold was entered in CASP14 under the team name ‘AlphaFold2’.
This system used a completely different model from our CASP13 entry31,
and demonstrated a considerable improvement over previous methods
in terms of providing routinely high accuracy29,30. Backbone predic-
tions with sub-Ångström root mean square deviation (Cα r.m.s.d.)
are now common, and side chains are increasingly accurate2. Good
results can often be achieved even for challenging proteins without a
template structure in the PDB, or with relatively few related sequences
to build a multiple sequence alignment (MSA)2.

We predicted structures for the UniProt human reference proteome
(one representative sequence per gene), with an upper length limit of
2,700 residues6. The final dataset covers 98.5% of human proteins with
a full chain prediction.

we expect to see high confidence on domains but low confidence on
linkers and unstructured regions (Extended Data Fig. 1). To this end,
AlphaFold produces a per-residue confidence metric called predicted
local distance difference test (pLDDT) on a scale from 0 to 100. pLDDT
estimates how well the prediction would agree with an experimental
structure based on the local distance difference test Cα (lDDT-Cα)35.
It has been shown to be well-calibrated (Fig. 1a, Extended Data Fig. 2
and Extended Data Table 1) and full details on how the pLDDT is pro-
duced are given in the supplementary information of the companion
AlphaFold paper2.
We consider a prediction highly accurate when—in addition to a
good backbone prediction—the side chains are frequently correctly
oriented. On this basis, pLDDT > 90 is taken as the high accuracy cut-off,
above which AlphaFold χ1 rotamers are 80% correct for a recent PDB

Of the human proteome, 35.7% of total residues fall within the high-
est accuracy band (corresponding to 38.6% of residues for which a
prediction was produced) (Fig. 1c). This is double the number of resi-
dues covered by an experimental structure. In total, 58.0% of residues
were predicted confidently (pLDDT > 70), indicating that we also add
substantial coverage for sequences without a good template in PDB
(with a sequence identity below 30%). At the per-protein level, 43.8% of
proteins have a confident prediction on at least three quarters of their
sequence, while 1,290 proteins contain a substantial region (more than
200 residues) with pLDDT ≥ 70 and no good template.

Membrane proteins, in particular, are generally underrepresented in the
PDB because they have historically been challenging experimental tar-
gets. This shows that AlphaFold is able to produce confident predictions
even for protein classes that are not abundant within its training set.

Many previous large-scale structure prediction efforts have focused on
domains—regions of the sequence that fold independently9–11,15. Here
we process full-length protein chains. There are several motivations for
this. Restricting the prediction to pre-identified domains risks miss-
ing structured regions that have yet to be annotated. It also discards
contextual information from the rest of the sequence, which might be
useful in cases in which two or more domains interact substantially.

Finally, the full chain approach lets the model attempt an inter-domain
packing prediction.
Inter-domain accuracy was assessed at CASP14, and AlphaFold out-
performed other methods41. However, the assessment was based on a
small target set. To further evaluate AlphaFold on long multi-domain
proteins, we compiled a test dataset of recent PDB chains that were
not in the training set of the model. Only chains with more than 800
resolved residues were included, and a template filter was applied
(Methods). Performance on this set was evaluated using the template
modelling score (TM-score42), which should better reflect global, as
opposed to per-domain, accuracy. The results were encouraging, with
70% of predictions having a TM-score > 0.7 (Fig. 2a).

In a systematic
analysis of recent PDB chains, we observed that AlphaFold has much
lower accuracy for regions in which the chain has a high percentage of
heterotypic, cross-chain contacts (Fig. 4d).
In summary, our current interpretation of regions in which AlphaFold
exhibits low pLDDT is that they have high likelihood of being unstruc-
tured in isolation.

The parts of the human proteome that are still without a confident
prediction represent directions for future research. Some proportion
of these will be genuine failures, in which a fixed structure exists but
the current version of AlphaFold does not predict it. In many other
cases, in which the sequence is unstructured in isolation, the problem
arguably falls outside the scope of single-chain structure prediction.

![fold3.png](https://github.com/LoqmanSamani/protein_sa/blob/systembiology/images/fold3.png)
![fold4.png](https://github.com/LoqmanSamani/protein_sa/blob/systembiology/images/fold4.png)
![fold5.png](https://github.com/LoqmanSamani/protein_sa/blob/systembiology/images/fold5.png)
![fold6.png](https://github.com/LoqmanSamani/protein_sa/blob/systembiology/images/fold6.png)






