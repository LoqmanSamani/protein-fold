## Highly accurate protein structure prediction with AlphaFold

![fold20](https://github.com/LoqmanSamani/protein_sa/blob/systembiology/images/fold20.png)


## 1 & 2 . Genetic and Structure Database Searches

The first step in the alphafold algorithm is **Data Pipeline**, which is the process of collecting input information for the rest of the alphafold process.
The initial input in this stage is an mmCIF (in training phase) or a fasta (in prediction phase) file which produces input features for the model.
Next the input file is parsed and relevant parameters is extracted (mmCIF: the sequence, atom coordinates, release date, name, and resolution. fasta: the sequence and its name).
with the information in hand two different database searches is conducted, in one case, genetic search, the database Big Fantastic Database (BFD), created using Uniprot/TrEMBL+Swissprot database,
is searched for similar sequences (in this case specific applications like JackHMMER and HHBlits are used) and the results, called multiple sequence alignments (MSAs), then are cleaned from duplicated sequences 
and the remaining sequences then stacked together to form a matrix of sequences. Each column of the MSA representation encodes an
individual residues of the input sequence while each row represent the sequences in which those residues appear.
in the next case, the information obtained from the last step (MSAs) are used (with an approach called HHSearch) to search in Protein Data Bank (PDB), template search,  to obtain the protein structures, which 
have a similar structure as the input MSA (in many cases the 3D structure and also the functionality of two protein in different organisms are similar but the primary structures may be different, so in this step the proteins with the same structure are selected). 
These structural information are used as pair representation (an initial representation of the 3D structure of the target protein), which contains the structural information in form a 3D matrix with the initial spatial information of the sequences.
The collected structural information (pair representations) are limited to 4 sequences (inference phase) or 20 sequences (training phase).

in training phase 75% of training data comes from genetic search (sequence information) approach and the rest 25% comes from template search (pair representations).

co-evolution: Suppose we have a protein where an amino acid with negative charge (say, glutamate) is near to an amino acid with positive charge (say, lysine), although they are both far away in the amino acid sequence. 
This Coulombic interaction stabilises the structure of the protein. Imagine now that the first amino acid mutates into a positively charged amino acid — in order to preserve this contact, 
the second amino acid will be under evolutionary pressure to mutate into a negatively charged amino acid, otherwise the resulting protein may not be able to fold.


## 3 . Evoformer (evolutionary transformer?) module

Evoformer is the main step in the alphafold algorithm, it is actually a transformer neural network which takes the MSAs and structural information from the last steps
as input and the output of this network is the prediction of protein structures as a graph inference problem in 3D space in which the edges of the
graph are defined by residues in proximity, it means each node in the graph corresponds to a residue (amino acid), and the edges between nodes represent the spatial relationships between those residues.
The objective of this part is to refine the representations for both the MSA and the pair interactions, but also to iteratively exchange information between them.
The Evoformer network consists of 48 (by default) blocks (iterations), each block has an MSA representation and a pair representation as its input and output and processes them within several layers.

First step in this network is embedding step (In AlphaFold 2, the embeddings are vanilla dense neural networks.), which means to convert the alphabetic information (A, C, G, T) stores in MSA and also in pair representation to numbers (probabilities), which are the input of the neural networks in general.
Input primary sequence and MSA features are embedded to MSA representation which is a matrix with the shape of (number of sequences * length of the sequence). The structural information are also embedded to form pair representation, a matrix with the shape of (length of the sequence * length of the sequence).


----
The objective of attention is to identify which parts of the input are more important for the objective of the neural network. In other words, to identify which parts of the input it should pay attention to.

GPT: Generative Pre-trained Transformer

There is a reason why transformers have not been widely implemented in many fields. The construction of the attention matrix leads to a quadratic memory cost. 
The Evoformer architecture uses not one, but two transformers (a “two-tower architecture”), with one clear communication channel between the two. Each head is specialised for the particular type of data it is looking at, either a multiple sequence alignment, or a representation of pairwise interactions between amino acids. They also incorporate the information of the contiguous representation, allowing for regular exchange of information and iterative refinement.

The MSA transformer computes attention over a very large matrix of protein symbols. To reduce what would otherwise be an impossible computational cost, the attention is “factorised” in “row-wise” and “column-wise” components. Namely, the network first computes attention in the horizontal direction, allowing the network to identify which pairs of amino acids are more related; and then in the vertical direction, determining which sequences are more informative.

The most important feature of AlphaFold 2’s MSA transformer is that the row-wise (horizontal) attention mechanism incorporates information from the “pair representation”. When computing attention, the network adds a bias term that is calculated directly from the current pair representation. This trick augments the attention mechanism and allows it to pinpoint interacting pairs of residues.

The other transformer head, the one that acts on the pair representation, works in a similar manner, although a lot of details differ, of course. The key feature of this network is that attention is arranged in terms of triangles of residues. The intuition here is to enforce the triangle inequality, one of the axioms of metric spaces. This is quite a clever idea since one of the classical problems of deep learning-based structure prediction was that distance distributions could not be embedded in three-dimensional space. It seems this trick fixes that and then some more.


----

## 4. The structure module

this part of the alphafold algorithm takes the refined MSA representation and pair representation, which were manipulated with Evoformer, and leverages them to construct a three-dimensional model of the structure.
The end result is a long list of Cartesian coordinates representing the position of each atom of the protein, including side chains.

## 5. Recicling

The last step of the alphafold algorithm is the Recycling step. After generating the final structure, output of th structure module, all the generated and refined information (i.e. MSA representation, pair representation and predicted structure) will be used as input for Evoformer blocks, This allows the model to refine its predictions.
This process is repeated 3 times (by default) and in each iteration the output is used as input of the model, so the prediction will be more precise and accurate.


1. MSA stack: 

The final Evoformer block provides a highly processed MSA representation {msi } and a pair representation {zij }, which contain information required for the structure module










































three main parts of the AlphaFold 2 system: 
[2.png]()

Transformer:  an “oracle” that can quickly identify which pieces of information are more informative.






1. number of residues in the input primary sequence by Nres (cropped during training). 
2. the number of templates used in the model by Ntempl.
3. the number of all available MSA sequences by Nall_seq.
4. the number of clusters after MSA clustering by Nclust.
5. the number of sequences processed in the MSA stack by Nseq (where Nseq = Nclust + Ntempl )
6. the number of unclustered MSA sequences by Nextra_seq.



-------------------------------------------------------------------------------------------------------
### I have a question about the coordinates of an atom in an  PDB file or in an mmCIF file:    which point will be used as (0, 0, 0) for this coordinates


In both PDB (Protein Data Bank) and mmCIF (macromolecular Crystallographic Information File) formats, the coordinate system is determined by the crystallographic or experimental methods used to determine the structure of the biomolecule. The origin (0, 0, 0) in the coordinate system is typically chosen to be a reference point within the crystal structure.

For PDB files, the origin is often set to the center of mass or another reference point in the crystallographic unit cell. The origin is chosen based on the crystallographic symmetry and the arrangement of biomolecules within the crystal lattice.

For mmCIF files, which are more structured and flexible than PDB files, the origin is also determined based on the crystallographic information. The mmCIF format provides a standardized way to represent crystallographic information, including unit cell parameters and space group details.

To determine the origin, you would need to refer to the experimental details in the header of the PDB or mmCIF file. Look for information such as crystallographic symmetry, unit cell parameters, and space group information. The exact method and choice of origin may vary depending on the specifics of the experiment.

Keep in mind that the coordinates in these files are given in angstroms or nanometers, depending on the convention used, and they represent the spatial positions of atoms within the crystal structure.


