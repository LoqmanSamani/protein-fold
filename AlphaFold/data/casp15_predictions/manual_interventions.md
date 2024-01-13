# Manual interventions (per target)

## T1123
We used 200 seeds on one of our standard models instead of the automated system
described above, though the prediction is not much different from the automated
prediction.

## T1125
Using updated genetics we found more sequences, but we also supplemented these
with 12 extra sequences found by running PsiBLAST against GenBank. The extra
sequences mostly covered a single domain, and the models still didn't produce a
convincing domain packing. Running 200 seeds each for multimer and single-chain
models, we found a slightly more plausible model among the single-chain
predictions.

## T1130 / T1131
These two targets yielded very few hits using standard MSA search methods. Based
on the species name, the targets appeared to be related to the proteins
described in this paper, which was available before the target closed:
https://academic.oup.com/gbe/article/14/6/evac069/6582783. Homologous sequences
were therefore sourced from the supplementary information files here:
https://janelia.figshare.com/articles/dataset/Multi-sequence_alignment_files/19560991,
with each target aligned to the MSA(s) containing a similar sequence. T1130 was
predicted using the new AlphaFold-Multimer models with the custom MSA. T1131 was
predicted using a custom MSA as well but just with the released AlphaFold model
(prediction was already much improved so didn’t try the new models).

## T1169
This target was run with the model trained with crop size 896 residues and 30
recycles instead of 640 residues and 20 recycles, which likely worked better
due to the very large single chain.

## H1111
This target was described as A9B9C9 but the full complex was too large for our
system, so we ran a few smaller stoichiometries, including A1B1C1, A2B2C2, and
A3B3C3. We then tried aligning the most confident predictions from each run with
the 9-mer PDB entry 7ALW, which was provided to predictors as a partial model in
the "Additional Information" section. By manually evaluating the quality of the
alignments and interfaces we chose the "middle" trimer of the A3B3C3 to copy 9
times. Each copy was aligned with one of the 9 chains in 7ALW to construct the
full A9B9C9 model. The assembled structure was then relaxed.

## H1106
This target was run with our automated system as an A1B1 complex, whose selected
conformation forms a slightly different conformation than as a subcomplex of
H1111. Interestingly, among our multiple seeds for H1106, there is a near even
mixture of two states where the other state looks like the state in H1111. We
are unsure if this indicates conformational change upon binding or just
prediction uncertainty.

## H1114
We were unable to run a full length model of this complex due to size. We first
ran A4B8 stoichiometry and observed a B2 dimer docking against each of the A
chains where the A chains formed a ring. Based on this, we ran a B2C2 complex
which produced B2 predictions consistent with the first run. We then aligned the
relaxed A4B8 model to four copies of the relaxed B2C2 prediction to make a full
complex model (taking the B2 chain from the B2C2 model to minimize
stereochemistry violations).

## H1137
Automated full complex predictions always curved back over in the tube region in
a way that looked implausible to us (though we might have submitted among 5
models if entering CASP).

Running the "tube" part of the protein in isolation produced a straight
prediction so we chose to assemble the prediction in three parts: the base
(transmembrane and GFP), the tube, and the post-tube region. We aligned separate
pre-relaxed predictions for each of these on overlapping regions of about 30
residues each, then ran the final Amber relax on the combined model to make the
final model.

---

In addition to the manual interventions above, the CASP15 baselines were run
with a different method to choose the number of MSA sequences for each chain
(up to 16384/num_chains per chain rather than fixed up to 2048). In general on a
large set of proteins, this new scheme doesn't seem to improve results and may
make them worse so we didn't include it in this release. We also had slightly
different HHBlits prefilters on our baseline run.

