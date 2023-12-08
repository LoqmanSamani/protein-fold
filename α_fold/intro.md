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

In particular, the latest version of
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


3. [Highly accurate protein structure prediction with AlphaFold](https://scholar.google.com/scholar?hl=en&as_sdt=0%2C5&q=Jumper%2C+J.+et+al.+Highly+accurate+protein+structure+prediction+with+AlphaFold.+Nature+https%3A%2F%2Fdoi.org%2F10.1038%2Fs41586-021-03819-2+%282021%29&btnG=)

Proteins are essential to life, and understanding their structure can facilitate a
mechanistic understanding of their function. Through an enormous experimental
effort1–4, the structures of around 100,000 unique proteins have been determined5, but
this represents a small fraction of the billions of known protein sequences6,7.
the structure prediction component of the ‘protein folding problem’8—has
been an important open research problem for more than 50 years9.

Here we provide the first computational method
that can regularly predict protein structures with atomic accuracy even in cases in which
no similar structure is known. We validated an entirely redesigned version of our neural
network-based model, AlphaFold, in the challenging 14th Critical Assessment of protein
Structure Prediction (CASP14)15, demonstrating accuracy competitive with
experimental structures in a majority of cases and greatly outperforming other
methods. Underpinning the latest version of AlphaFold is a novel machine learning
approach that incorporates physical and biological knowledge about protein structure,
leveraging multi-sequence alignments, into the design of the deep learning algorithm.

The development of computational methods to predict
three-dimensional (3D) protein structures from the protein sequence
has proceeded along two complementary paths that focus on either the
physical interactions or the evolutionary history. The physical interac-
tion programme heavily integrates our understanding of molecular
driving forces into either thermodynamic or kinetic simulation of pro-
tein physics16 or statistical approximations thereof17. Although theoreti-
cally very appealing, this approach has proved highly challenging for
even moderate-sized proteins due to the computational intractability
of molecular simulation, the context dependence of protein stability
and the difficulty of producing sufficiently accurate models of protein
physics. The evolutionary programme has provided an alternative in
recent years, in which the constraints on protein structure are derived
from bioinformatics analysis of the evolutionary history of proteins,
homology to solved structures18,19 and pairwise evolutionary correla-
tions20–24. This bioinformatics approach has benefited greatly from
the steady growth of experimental protein structures deposited in
the Protein Data Bank (PDB)5,

In this study, we develop the first, to our knowledge, computational
approach capable of predicting protein structures to near experimental
accuracy in a majority of cases.

The CASP assessment is
carried out biennially using recently solved structures that have not
been deposited in the PDB or publicly disclosed so that it is a blind test
for the participating methods, and has long served as the gold-standard
assessment for the accuracy of structure prediction25,26.
AlphaFold structures had a median backbone
accuracy of 0.96 Å r.m.s.d.95 (Cα root-mean-square deviation at 95%
residue coverage) (95% confidence interval = 0.85–1.16 Å) whereas
the next best performing method had a median backbone accuracy
of 2.8 Å r.m.s.d.95 (95% confidence interval = 2.7–4.0 Å) (measured on
CASP domains;

As a comparison point for this accuracy,
the width of a carbon atom is approximately 1.4 Å.

. The all-atom accuracy of Alpha-
Fold was 1.5 Å r.m.s.d.95 (95% confidence interval = 1.2–1.6 Å) compared
with the 3.5 Å r.m.s.d.95 (95% confidence interval = 3.1–4.2 Å) of the best
alternative method. Our methods are scalable to very long proteins with
accurate domains and domain-packing (see Fig. 1d for the prediction
of a 2,180-residue protein with no structural homologues). F

We demonstrate in Fig. 2a that the high accuracy that AlphaFold dem-
onstrated in CASP14 extends to a large sample of recently released PDB
structures; in this dataset, all structures were deposited in the PDB after
our training data cut-off and are analysed as full chains.

The network comprises two main stages. First, the trunk of the net-
work processes the inputs through repeated layers of a novel neural
network block that we term Evoformer to produce an Nseq × Nres array
(Nseq, number of sequences; Nres, number of residues) that represents
a processed MSA and an Nres × Nres array that represents residue pairs.
The MSA representation is initialized with the raw MSA (although
see Supplementary Methods 1.2.7 for details of handling very deep
MSAs). The Evoformer blocks contain a number of attention-based
and non-attention-based components. We show evidence in ‘Interpret-
ing the neural network’ that a concrete structural hypothesis arises
early within the Evoformer blocks and is continuously refined. The key
innovations in the Evoformer block are new mechanisms to exchange
information within the MSA and pair representations that enable direct
reasoning about the spatial and evolutionary relationships.
The trunk of the network is followed by the structure module that
introduces an explicit 3D structure in the form of a rotation and transla-
tion for each residue of the protein (global rigid body frames). These
representations are initialized in a trivial state with all rotations set to
the identity and all positions set to the origin, but rapidly develop and
refine a highly accurate protein structure with precise atomic details.
Key innovations in this section of the network include breaking the
chain structure to allow simultaneous local refinement of all parts of
the structure, a novel equivariant transformer to allow the network to
implicitly reason about the unrepresented side-chain atoms and a loss
term that places substantial weight on the orientational correctness
of the residues. Both within the structure module and throughout
the whole network, we reinforce the notion of iterative refinement
by repeatedly applying the final loss to outputs and then feeding the
outputs recursively into the same modules. The iterative refinement
using the whole network (which we term ‘recycling’ and is related to
approaches in computer vision28,29) contributes markedly to accuracy
with minor extra training time (see Supplementary Methods 1.8 for
details).

Predictions of side-chain χ angles as well as the final, per-residue
accuracy of the structure (pLDDT) are computed with small per-residue
networks on the final activations at the end of the network. The estimate
of the TM-score (pTM) is obtained from a pairwise error prediction that
is computed as a linear projection from the final pair representation. The
final loss (which we term the frame-aligned point error (FAPE) (Fig. 3f))
compares the predicted atom positions to the true positions under
many different alignments. For each alignment, defined by aligning
the predicted frame (Rk, tk) to the corresponding true frame, we com-
pute the distance of all predicted atom positions xi from the true atom
positions. The resulting Nframes × Natoms distances are penalized with a
clamped L1 loss. This creates a strong bias for atoms to be correct relative
to the local frame of each residue and hence correct with respect to its
side-chain interactions, as well as providing the main source of chirality
for AlphaFold (Supplementary Methods 1.9.3 and Supplementary Fig. 9).

The AlphaFold architecture is able to train to high accuracy using only
supervised learning on PDB data, but we are able to enhance accuracy
(Fig. 4a) using an approach similar to noisy student self-distillation35.
In this procedure, we use a trained network to predict the structure of
around 350,000 diverse sequences from Uniclust3036 and make a new
dataset of predicted structures filtered to a high-confidence subset. We
then train the same architecture again from scratch using a mixture of
PDB data and this new dataset of predicted structures as the training
data, in which the various training data augmentations such as crop-
ping and MSA subsampling make it challenging for the network to
recapitulate the previously predicted structures. This self-distillation
procedure makes effective use of the unlabelled sequence data and
considerably improves the accuracy of the resulting network.
Additionally, we randomly mask out or mutate individual residues
within the MSA and have a Bidirectional Encoder Representations from
Transformers (BERT)-style37 objective to predict the masked elements of
the MSA sequences. This objective encourages the network to learn to
interpret phylogenetic and covariation relationships without hardcoding
a particular correlation statistic into the features. The BERT objective is
trained jointly with the normal PDB structure loss on the same training
examples and is not pre-trained, in contrast to recent independent work38.

To understand how AlphaFold predicts protein structure, we trained
a separate structure module for each of the 48 Evoformer blocks in
the network while keeping all parameters of the main network fro-
zen (Supplementary Methods 1.14). Including our recycling stages,
this provides a trajectory of 192 intermediate structures—one per full
Evoformer block—in which each intermediate represents the belief of
the network of the most likely structure at that block. The resulting
trajectories are surprisingly smooth after the first few blocks, show-
ing that AlphaFold makes constant incremental improvements to the
structure until it can no longer improve (see Fig. 4b for a trajectory of
accuracy). These trajectories also illustrate the role of network depth.
For very challenging proteins such as ORF8 of SARS-CoV-2 (T1064),
the network searches and rearranges secondary structure elements
for many layers before settling on a good structure. For other proteins
such as LmrP (T1024), the network finds the final structure within the
first few layers.

Although AlphaFold has a high accuracy across the vast majority of
deposited PDB structures, we note that there are still factors that affect
accuracy or limit the applicability of the model. The model uses MSAs
and the accuracy decreases substantially when the median alignment
depth is less than around 30 sequences (see Fig. 5a for details). We
observe a threshold effect where improvements in MSA depth over
around 100 sequences lead to small gains. We hypothesize that the MSA
information is needed to coarsely find the correct structure within the
early stages of the network, but refinement of that prediction into a
high-accuracy model does not depend crucially on the MSA information.

These approaches effectively leverage the rapid improvement in com-
puter vision systems46 by treating the problem of protein structure
prediction as converting an ‘image’ of evolutionary couplings22–24 to an
‘image’ of the protein distance matrix and then integrating the distance
predictions into a heuristic system that produces the final 3D coordinate
prediction.

![fold7.png](https://github.com/LoqmanSamani/protein_sa/blob/systembiology/images/fold7.png)
![fold8.png](https://github.com/LoqmanSamani/protein_sa/blob/systembiology/images/fold8.png)
![fold9.png](https://github.com/LoqmanSamani/protein_sa/blob/systembiology/images/fold9.png)
![fold10.png](https://github.com/LoqmanSamani/protein_sa/blob/systembiology/images/fold7.png)



4. [Deep learning techniques have significantly impacted protein structure prediction and protein design](https://scholar.google.com/scholar?hl=en&as_sdt=0%2C5&q=Deep+learning+techniques+have+significantly+impacted+protein+structure+prediction+and+protein+design+Robin+Pearce1+and+Yang+Zhang&btnG=)

Given the cost, both financially and
timewise, associated with experimentally determining a
protein’s structure and function, extensive effort has been
made to develop computational methods capable of
modeling the structures of natural protein sequences
and/or designing new sequences with novel structures
and functions beyond proteins observed in nature.

The goal of protein structure prediction is to use compu-
tational methods to determine the spatial location of
every atom in a protein molecule starting from its amino
acid sequence. Depending on whether a template struc-
ture is used, protein structure prediction approaches can
be generally categorized as either template-based model-
ing (TBM) or template-free modeling (FM) methods.
While TBM constructs models by copying and refining
structural frameworks of other related proteins, called
templates, identified from the PDB, FM aims to
predict protein structures without using global template
structures. FM methods have also been referred to as ab
initio or de novo modeling approaches. A general pipeline
that illustrates the key steps involved in traditional TBM
and FM methods is depicted in Figure 1.

There are four key steps involved in TBM methods: (1)
identification of experimentally solved proteins (tem-
plates) structurally related to the protein to be modeled
(query), (2) alignment of the query and the template
proteins, (3) construction of the initial structure frame-
works by copying the aligned regions of the template
structure, and (4) construction of the unaligned regions
and refinement of the global structure.

Depending on the evolutionary distance between the
query and template, TBM has been historically divided
into comparative modeling (CM), which is designed for
targets with close homologous templates where the tem-
plates can typically be identified by sequence-based
alignment, and threading, which is designed for detecting
more distantly homologous templates by combining
sequence profiles and/or Hidden Markov Model
(HMM) alignment with local structure feature prediction
[2,3].

One of the earliest
sequence-based contact prediction methods used corre-
lated mutations observed in multiple sequence align-
ments (MSAs) to predict inter-residue contact maps
[21]. The hypothesis behind the approach was that if
mutations that occur at two positions in an MSA are
correlated, these positions are more likely to form a
contact in 3D space. This is because there is evolutionary
pressure to conserve the structures of proteins and a
mutation at one position may be rescued by a correspond-
ing mutation at a nearby residue. The accuracy of co-
evolution-based contact map prediction remained low for
many years due to the inability to distinguish between
direct and indirect interactions, where indirect interac-
tions occur when residues appear to co-evolve but do not
actually form contacts. For example, if Residues A and B
are both in contact with Residue C, A and B often appear
as if they co-evolve even when there is no physical
contact between them. There is evidence showing that
such co-evolution may have a functional cause [22] rather
than a structural one, which resulted in the failure of
structure-based contact derivation.

In addition to the high accuracy of
model training enabled by multi-layer neural networks
[1], another important advantage of deep learning is its
ability to predict multiple structural features, including
contacts, distances, inter-residue torsion angles and
hydrogen bonds. The combination of these structural
features with the classical folding simulation methods
has significantly improved the modeling accuracy of
protein structure prediction, especially for FM protein
targets which lack homologous templates.

The success of this approach can be partially attributed
to the ability of deep learning to simultaneously consider
the global set of pair-wise interactions instead of consider-
ing only one interaction at a time, thereby leading to more
accurate discrimination between direct and indirect con-
tacts .

A similar residual neural network was later extended to
predict the probability that the distance between two
residues falls within a given distance range instead of
predicting a binary contact map [32]. The power of
distance map-guided folding was convincingly demon-
strated by AlphaFold in the CASP13 experiment, in
which the program utilized an ultra-deep neural network
composed of 220 residual blocks to predict distance maps
for a query sequence [19]. The distance maps were then
used to guide their fragment assembly and gradient
descent-based folding simulations for full-length struc-
ture construction. AlphaFold also used a unique fragment
generation strategy where they leveraged deep learning
to produce short structural fragments de novo. To accom-
plish this, they trained a generative network to create
fragments based on prediction of the torsional angles for a
selected region of a protein. This approach allows for the
generation of fragments conditioned on the input features
and eliminates the need to identify near native fragments
from a library of existing fragment structures.

The most exciting progress in the history of protein
structure prediction was recently brought about by Alpha-
Fold2, the second iteration of AlphaFold developed by
the Google DeepMind team [35], which achieved an
unprecedented modeling accuracy in the CASP14 exper-
iment. Out of the 89 domains with experimentally
released structures, AlphaFold2 generated first-rank
models with TM-scores >0.5 for 88 domains, where
59 of them had TM-scores >0.914. Here, TM-score is
a sequence length-independent metric that measures
protein structural similarity and takes a value in the range
of (0, 1) [36], where PDB statistics show that a TM-score >
0.5 indicates that two structures share approximately
the same SCOP/CATH fold [37]. Moreover, we collected
a set of 112 single-domain proteins whose structures were
solved by both NMR and X-ray crystallography and had
sequence identities >95% and alignment gaps <10 AA,
where we found the average TM-score was 0.807 
0.107 between the NMR and X-ray structures. Thus,
AlphaFold2 could fold nearly all individual domains in
CASP14, with around 2/3 (=59/89) of the cases having
accuracy comparable to low-to-medium resolution exper-
imental models if we use a cutoff TM-score of 0.914
(=0.807 + 0.107).

These
data suggest that AlphaFold2 nearly solved the problem
of single-domain protein structure prediction, at least at
the fold level.

AlphaFold2 outperformed the second-best group by
a large margin with the average TM-score differing by
23% (0.903 versus 0.732).

Compared to the first iteration of AlphaFold in CASP13,
which was driven by convolutional neural network-
based distance map prediction, one of the major new
developments of AlphaFold2 is the attention-based
neural network architecture that attends arbitrarily over
the full MSA, which allows the system to select relevant
sequences from the MSAs and extract richer input
information. Moreover, instead of using gradient
descent optimization to construct models based on
the predicted distance restraints, as AlphaFold did in
CASP13, AlphaFold2 utilizes a full end-to-end training
system from sequence to structure models using itera-
tive structural refinement based on local structural
error estimation. As part of this, the system replaces
traditional folding simulations with a structure module
composed of 3D equivariant transformer neural net-
works, which treat each amino acid as a gas of 3D rigid
bodies and directly builds the protein backbone and
side-chains. All these advantages, together with the
extensive computing resources that are beyond what
are accessible to most of the academic research labora-
tories, contribute to the significant improvement of the
state-of-the-art of deep learning-based protein structure
prediction

The prediction of protein structures from amino acid
sequences alone has remained an outstanding problem
in structural biology since Anfisen first demonstrated that
the information encoded in a protein sequence determines
its structure more than 60 years ago. Now more than ever,
there is an urgent need to develop high accuracy
protein structure prediction methods, as advancements
in high-throughput sequencing technology have greatly
exacerbated the gap between the number of known
sequences and the number of experimentally determined
protein structures.

![fold11.png](https://github.com/LoqmanSamani/protein_sa/blob/systembiology/images/fold11.png)


5. [Improved protein structure prediction using potentials from deep learning](https://scholar.google.com/scholar?hl=en&as_sdt=0%2C5&q=Improved+protein+structure+prediction+using+potentials+from+deep+learning+Andrew+W.+Senior1%2C4*%2C+Richard+Evans1%2C4%2C+John+Jumper1%2C4%2C+James+Kirkpatrick1%2C4%2C+Laurent+Sifre1%2C4%2C+Tim+Green1%2C+Chongli+Qin1%2C+Augustin+%C5%BD%C3%ADdek1%2C+Alexander+W.+R.+Nelson1%2C+Alex+Bridgland1%2C+Hugo+Penedones1%2C+Stig+Petersen1%2C+Karen+Simonyan1%2C+Steve+Crossan1%2C+Pushmeet+Kohli1%2C+David+T.+Jones2%2C3%2C+David+Silver1%2C+Koray+Kavukcuoglu1+%26+Demis+Hassabis1&btnG=)


Protein structure prediction can be used to determine the three-dimensional shape of
a protein from its amino acid sequence1. This problem is of fundamental importance
as the structure of a protein largely determines its function2; however, protein
structures can be difficult to determine experimentally.
It is possible to infer which
amino acid residues are in contact by analysing covariation in homologous
sequences, which aids in the prediction of protein structures3.
We find that the resulting potential can be
optimized by a simple gradient descent algorithm to generate structures without
complex sampling procedures.
In the recent Critical
Assessment of Protein Structure Prediction5 (CASP13)—a blind assessment of the state
of the field—AlphaFold created high-accuracy structures (with template modelling
(TM) scores6 of 0.7 or higher) for 24 out of 43 free modelling domains, whereas the
next best method, which used sampling and contact information, achieved such
accuracy for only 14 out of 43 domains.

As the function of
a protein is dependent on its structure, understanding protein struc-
tures has been a grand challenge in biology for decades. Although
several experimental structure determination techniques have been
developed and improved in accuracy, they remain difficult and time-
consuming2.

CASP5 is a biennial blind protein structure prediction assessment
run by the structure prediction community to benchmark progress in
accuracy. In 2018, AlphaFold joined 97 groups from around the world in
entering CASP138. Each group submitted up to 5 structure predictions
for each of 84 protein sequences for which experimentally determined
structures were sequestered. Assessors divided the proteins into 104
domains for scoring and classified each as being amenable to template-
based modelling (TBM, in which a protein with a similar sequence has
a known structure, and that homologous structure is modified in
accordance with the sequence differences) or requiring free model-
ling (FM, in cases in which no homologous structure is available), with
an intermediate (FM/TBM) category. Figure 1a shows that AlphaFold
predicts more FM domains with high accuracy than any other system,
particularly in the 0.6–0.7 TM-score range. The TM score—ranging
between 0 and 1—measures the degree of match of the overall (back-
bone) shape of a proposed structure to a native structure. The assessors
ranked the 98 participating groups by the summed, capped z-scores of
the structures, separated according to category. AlphaFold achieved
a summed z-score of 52.8 in the FM category (best-of-five) compared
with 36.6 for the next closest group (322). Combining FM and TBM/FM
categories, AlphaFold scored 68.3 compared with 48.2.

Despite using only FM techniques and not using templates, AlphaFold
also scored well in the TBM category according to the assessors’ for-
mula 0-capped z-score, ranking fourth for the top-one model or first
for the best-of-five models. Much of the accuracy of AlphaFold is due
to the accuracy of the distance predictions, which is evident from the
high precision of the corresponding contact predictions.

We show
that it is possible to construct a learned, protein-specific potential
by training a neural network (Fig. 2b) to make accurate predictions
about the structure of the protein given its sequence, and to predict
the structure itself accurately by minimizing the potential by gradient
descent (Fig. 2c). The neural network predictions include backbone
torsion angles and pairwise distances between residues. Distance
predictions provide more specific information about the structure
than contact predictions and provide a richer training signal for the
neural network. By jointly predicting many distances, the network
can propagate distance information that respects covariation, local
structure and residue identities of nearby residues. The predicted
probability distributions can be combined to form a simple, principled
protein-specific potential. We show that with gradient descent, it is
simple to find a set of torsion angles that minimizes this protein-specific
potential using only limited sampling. We also show that whole chains
can be optimized simultaneously, avoiding the need to segment long
proteins into hypothesized domains that are modelled independently
as is common practice.

The central component of AlphaFold is a convolutional neural
network that is trained on PDB structures to predict the distances
dij between the Cβ atoms of pairs, ij, of residues of a protein. On the
basis of a representation of the amino acid sequence, S, of a protein
and features derived from the MSA(S) of that sequence, the network,
which is similar in structure to those used for image-recognition tasks29,
predicts a discrete probability distribution P(dij|S, MSA(S)) for every
ij pair in any 64 × 64 region of the L × L distance matrix, as shown in
Fig. 2b. The full set of distance distribution predictions constructed
by combining such predictions that covers the entire distance map is
termed a distogram (from distance histogram). Example distogram
predictions for one CASP protein, T0955, are shown in Fig. 3c, d. 

The high accuracy of the distance predictions
and consequently the contact predictions (Fig. 1c) comes from a com-
bination of factors in the design of the neural network and its training,
data augmentation, feature representation, auxiliary losses, cropping
and data curation.

We parameterized protein structures by the backbone
torsion angles (φ, ψ) of all residues and build a differentiable model of
protein geometry x = G(φ, ψ) to compute the Cβ coordinates, xi for all
residues i and thus the inter-residue distances, dij = ||xi − xj||, for each
structure, and express Vdistance as a function of φ and ψ. For a protein with
L residues, this potential accumulates L2 terms from marginal distribu-
tion predictions. To correct for the overrepresentation of the prior, we
subtract a reference distribution30 from the distance potential in the log
domain. The reference distribution models the distance distributions
P(dij|length) independent of the protein sequence and is computed
by training a small version of the distance prediction neural network
on the same structures, without sequence or MSA input features.
A separate output head of the contact prediction network is trained to
predict discrete probability distributions of backbone torsion angles
P(φi,ψi|S, MSA(S)). After fitting a von Mises distribution, this is used to
add a smooth torsion modelling term, Vtorsion, to the potential. Finally,
to prevent steric clashes, we add the Vscore2_smooth score of Rosetta9 to the
potential, as this incorporates a van der Waals term. We used multipli-
cative weights for each of the three terms in the potential; however, no
combination of weights noticeably outperformed equal weighting.
As all of the terms in the combined potential Vtotal(φ, ψ) are
differentiable functions of (φ, ψ), it can be optimized with respect to
these variables by gradient descent. Here we use L-BFGS31. Structures
are initialized by sampling torsion values from P(φi, ψi|S, MSA(S)).
Figure 2c illustrates a single gradient descent trajectory that minimizes
the potential, showing how this greedy optimization process leads to
increasing accuracy and large-scale conformation changes. The sec-
ondary structure is partly set by the initialization from the predicted
torsion angle distributions. The overall accuracy (TM score) improves
quickly and after a few hundred steps of gradient descent the accuracy
of the structure has converged to a local optimum of the potential.

We repeated the optimization from sampled initializations,
leading to a pool of low-potential structures from which further struc-
ture initializations are sampled, with added backbone torsion noise
(‘noisy restarts’), leading to more structures to be added to the pool.
After only a few hundred cycles, the optimization converges and the
lowest potential structure is chosen as the best candidate structure.
Figure 2e shows the progress in the accuracy of the best-scoring struc-
tures over multiple restarts of the gradient descent process, show-
ing that after a few iterations the optimization has converged. Noisy
restarts enable structures with a slightly higher TM score to be found
than when continuing to sample from the predicted torsion distribu-
tions (average of 0.641 versus 0.636 on our test set, shown in Extended
Data Fig. 4).

The following tools and dataset versions were used for the CASP sys-
tem and for subsequent experiments: PDB 15 March 2018; CATH 16
March 2018; HHblits based on v.3.0-beta.3 (three iterations, E = 1 × 10−3);
HHpred web server; Uniclust30 2017-10; PSI-BLAST v.2.6.0 nr dataset
(as of 15 December 2017) (three iterations, E = 1 × 10−3); SST web server
(March 2019); BioPython v.1.65; Rosetta v.3.5; PyMol 2.2.0 for structure
visualization; TM-align 20160521.

Our models are trained on structures extracted from the PDB13.
We extract non-redundant domains by utilizing the CATH34 35%
sequence similarity cluster representatives. This generated 31,247
domains, which were split into train and test sets (29,427 and 1,820
proteins, respectively).

For each training sequence, we searched for and aligned to the train-
ing sequence similar protein sequences in the Uniclust3035 dataset
with HHblits36 and used the returned MSA to generate profile features
with the position-specific substitution probabilities for each residue
as well as covariation features—the parameters of a regularized pseu-
dolikelihood-trained Potts model similar to CCMpred16. CCMPred uses
the Frobenius norm of the parameters, but we feed both this norm
(1 feature) and the raw parameters (484 features) into the network for
each residue pair ij. In addition, we provide the network with features
that explicitly represent gaps and deletions in the MSA. To make the
network better able to make predictions for shallow MSAs, and as a
form of data augmentation, we take a sample of half the sequences
from the the HHblits MSA before computing the MSA-based features.
Our training set contains 10 such samples for each domain. We extract
additional profile features using PSI-BLAST37.
The distance prediction neural network was trained with the follow-
ing input features (with the number of features indicated in brackets).
• Number of HHblits alignments (scalar).
• Sequence-length features: 1-hot amino acid type (21 features);
profiles: PSI-BLAST (21 features), HHblits profile (22 features),
non-gapped profile (21 features), HHblits bias, HMM profile (30
features), Potts model bias (22 features); deletion probability (1 fea-
ture); residue index (integer index of residue number, consecutive
except for multi-segment domains, encoded as 5 least-significant
bits and a scalar).
• Sequence-length-squared features: Potts model parameters
(484 features, fitted with 500 iterations of gradient descent using
Nesterov momentum 0.99,

The inter-residue distances are predicted by
a deep neural network. The architecture is a deep two-dimensional
dilated convolutional residual network.
The network is trained with stochastic gradient descent using a
cross-entropy loss. The target is a quantification of the distance
between the Cβ atoms of the residues (or Cα for glycine). We divide
the range 2–22 Å into 64 equal bins. The input to the network consists
of a two-dimensional array of features in which each i,j feature is the
concatenation of the one-dimensional features for both i and j as well
as the two-dimensional features for i,j.
Individual training runs were cross-validated with early stopping
using 27 CASP11 FM domains as a validation set. Models were selected
by cross-validation on 27 CASP12 FM domains.


![fold12.png]()
![fold13.png]()
![fold14.png]()
![fold15.png]()
![fold16.png]()




















