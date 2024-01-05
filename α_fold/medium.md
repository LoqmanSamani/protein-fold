## Highly accurate protein structure prediction with AlphaFold

![alphafold](https://github.com/LoqmanSamani/protein_sa/blob/systembiology/%CE%B1_fold/images/alphafold.png)


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



![evoformer]()

## 3 . Evoformer (evolutionary transformer?) module

Evoformer is the main step in the alphafold algorithm, it is actually a transformer neural network which takes the MSAs and structural information from the last steps
as input and the output of this network is the prediction of protein structures as a graph inference problem in 3D space in which the edges of the
graph are defined by residues in proximity, it means each node in the graph corresponds to a residue (amino acid), and the edges between nodes represent the spatial relationships between those residues.
The objective of this part is to refine the representations for both the MSA and the pair interactions, but also to iteratively exchange information between them.
The Evoformer network consists of 48 (by default) blocks (iterations), each block has an MSA representation and a pair representation as its input and output and processes them within several layers.

First step in this network is embedding step (In AlphaFold 2, the embeddings are vanilla dense neural networks.), which means to convert the alphabetic information (A, C, G, T) stores in MSA and also in pair representation to numbers (probabilities), which are the input of the neural networks in general.
Input primary sequence and MSA features are embedded to MSA representation which is a matrix with the shape of (number of sequences * length of the sequence). The structural information are also embedded to form pair representation, a matrix with the shape of (length of the sequence * length of the sequence).

- axial self-attention in the MSA stack
- triangular multiplicative updates and triangular self-attention in the pair stack
- outer product mean and attention biasing to allow communication between the stacks
- Each layer output is added via a residual connection to the current representations.
- Some layer outputs are passed through Dropout before they are added.(overfitting is a serious problem in such networks. Large networks are also slow to use, making it difficult to deal with overfitting by combining the predictions of many different large neural nets at test time. Dropout is a technique for addressing this problem. The key idea is to randomly drop units (along with their connections) from the neural network during training. This prevents units from co-adapting too much)

- input of Evoformer:  MSA representation {Msi}
- output: pair representation {Zij}


The final Evoformer block provides a highly processed MSA representation {msi } and a pair representation {zij }, which contain information required for the structure module.



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


----
AlphaFold 2 processes the sequence search data to generate two “representations”: a representation of the multiple sequence alignment (MSA), which captures sequence variation; and a representation of the “pairs of residues”, which captures which residues are likely to interact with each other.
The structure module considers the protein as a “residue gas”. Every amino acid is modelled as a triangle, representing the three atoms of the backbone. These triangles float around in space, and are moved by the network to form the structure.
These transformations are parametrised as “affine matrices”, which are a mathematical way to represent translations and rotations in a single 4×4 matrix:
At the beginning of the structure module, all of the residues are placed at the origin of coordinates. At every step of the iterative process, AlphaFold 2 produces a set of affine matrices that displace and rotate the residues in space. This representation does not reflect any physical or geometrical assumptions, and as a result the network has a tendency to generate structural violations
The secret sauce of the structure module is a new flavour of attention devised specifically for working with three-dimensional structures — DeepMind calls it “Invariant Point Attention” (IPA).
if the model knows that any possible rotation of translation of the data will lead to the same answer, it will need a lot less data to pull it away from wrong models and will therefore be able to learn much more.
 it is based on a very simple mathematical fact: that the L2-norm of a vector is invariant with respect to translations and rotations. 
One of said details is the loss function used by AlphaFold 2. The DeepMind team introduced a specific structural loss which they called FAPE (Frame Aligned Point Error), which we could understand as a clever version of the more commonly used root-mean-squared deviation (RMSD) of atomic positions.
It seems this function was not enough. The final AlphaFold 2 loss function is in fact a weighted sum of multiple “auxiliary losses”, which do not necessarily relate to the performance of the model, but that provide additional information. For example, the DeepMind team computed the loss not only of the final output, but also of each of the three iterations of the structure module. It also includes a “distogram loss”, where the predicted structure is used to generate a 2D matrix of distances that is also compared to the ground truth.
That is not the only sleight-of-hand. One of the cleverer tricks the team pulled off is self-distillation. In this approach, they took a model trained exclusively on the PDB, and predicted the structures of ~300k diverse protein sequences obtained from Uniclust. They then retrained the full model, incorporating a small random sample of these structures at every training cycle. They claim this allows the model to leverage the large amount of unlabelled data available in protein sequence repositories. 

----

## 5. Recycling

The last step of the alphafold algorithm is the Recycling step. After generating the final structure, output of th structure module, all the generated and refined information (i.e. MSA representation, pair representation and predicted structure) will be used as input for Evoformer blocks, This allows the model to refine its predictions.
This process is repeated 3 times (by default) and in each iteration the output is used as input of the model, so the prediction will be more precise and accurate.


1. MSA stack: 

The final Evoformer block provides a highly processed MSA representation {msi } and a pair representation {zij }, which contain information required for the structure module




 I realised what was the secret sauce:
It is their access to compute resources, and their top-notch know-how, that turned them into the successful neural network it became.

This leaves one to wonder exactly how much computational power went into this project. In their original announcement, DeepMind claimed that they used “128 TPUv3 cores or roughly equivalent to ~100-200 GPUs”. Although this amount of compute seems beyond the wildest dreams of most academic researchers, my friends at Google tell me it is not uncommon for individual Googlers to access similar resources on a regular basis. How many times more computational power did this team have, in comparison with all the other CASP14 contenders combined?
There won’t be an Android app anytime soon, but anyone should be able to run it, provided they are willing to invest in a powerful compute server (or the cloud).

if you have access to standard computing servers, then you should be able to run it without much trouble. I have tested it on a Quadro RTX 6000 GPU, (~20 GB of dedicated memory) and on a CPU server with ~300 GB of RAM, and in both cases I was able to obtain structure predictions. In particular, I have run a bunch of proteins of up to 600 amino acids, and I seem to have been able to produce an answer in every case.

the DeepMind team suggests that “accuracy drops substantially when the mean alignment depth is less than ~30 sequences”.

The most obvious direct application of this project is structure-based drug discovery. Until now, the availability of structures was a prime requirement: most people would not even consider starting a project without at least a crystal structure of the unbound protein. However, the availability of high-accuracy predictions, as well as predicted “error bars” will probably encourage the pharmaceutical industry to increasingly use AlphaFold 2’s models for development. And so, we may soon have inhibitors of many drug targets that have remained hitherto unexplored.


---
Transformer: a model that uses attention to boost the speed with which these models can be trained.
structure:
1. In a machine translation application, it would take a sentence in one language, and output its translation in another.
2. the transformer contains an encoding component, a decoding component, and connections between them.
3. The encoding component is a stack of encoders (6, 7 or ...). The decoding component is a stack of decoders of the same number.
4. The encoders are all identical in structure (yet they do not share weights). Each one is broken down into two sub-layers: a) self-attention b) feed forward neural network
      - self-attention layer helps the encoder look at other words in the input sentence as it encodes a specific word, the outputs of the self-attention layer are fed to a feed-forward neural network.
      - The decoder has both those layers, but between them is an attention layer that helps the decoder focus on relevant parts of the input sentence

how it works:
1. As is the case in NLP applications in general, we begin by turning each input word into a vector using an embedding algorithm.
2. The embedding only happens in the bottom-most encoder. The abstraction that is common to all the encoders is that they receive a list of vectors each of the size 512 – In the bottom encoder that would be the word embeddings, but in other encoders, it would be the output of the encoder that’s directly below. The size of this list is hyperparameter we can set – basically it would be the length of the longest sentence in our training dataset.
3. self-attention: ”The animal didn't cross the street because it was too tired”  What does “it” in this sentence refer to? Is it referring to the street or to the animal? It’s a simple question to a human, but not as simple to an algorithm. When the model is processing the word “it”, self-attention allows it to associate “it” with “animal”. As the model processes each word (each position in the input sequence), self attention allows it to look at other positions in the input sequence for clues that can help lead to a better encoding for this word.
   - The first step in calculating self-attention is to create three vectors from each of the encoder’s input vectors (in this case, the embedding of each word). So for each word, we create a Query vector, a Key vector, and a Value vector. These vectors are created by multiplying the embedding by three matrices that we trained during the training process.
   - Their dimensionality is 64, while the embedding and encoder input/output vectors have dimensionality of 512. They don’t HAVE to be smaller, this is an architecture choice to make the computation of multiheaded attention (mostly) constant.
   - The second step in calculating self-attention is to calculate a score. Say we’re calculating the self-attention for the first word in this example, “Thinking”. We need to score each word of the input sentence against this word. The score determines how much focus to place on other parts of the input sentence as we encode a word at a certain position.
   - The score is calculated by taking the dot product of the query vector with the key vector of the respective word we’re scoring. So if we’re processing the self-attention for the word in position #1, the first score would be the dot product of q1 and k1. The second score would be the dot product of q1 and k2.
   - The third and fourth steps are to divide the scores by 8 (the square root of the dimension of the key vectors used in the paper – 64. This leads to having more stable gradients. There could be other possible values here, but this is the default), then pass the result through a softmax operation. Softmax normalizes the scores so they’re all positive and add up to 1.
   - The fifth step is to multiply each value vector by the softmax score (in preparation to sum them up). The intuition here is to keep intact the values of the word(s) we want to focus on, and drown-out irrelevant words (by multiplying them by tiny numbers like 0.001, for example).
   - The sixth step is to sum up the weighted value vectors. This produces the output of the self-attention layer at this position (for the first word).
   - The paper further refined the self-attention layer by adding a mechanism called “multi-headed” attention (It expands the model’s ability to focus on different positions. It gives the attention layer multiple “representation subspaces”. As we’ll see next, with multi-headed attention we have not only one, but multiple sets of Query/Key/Value weight matrices (the Transformer uses eight attention heads, so we end up with eight sets for each encoder/decoder). Each of these sets is randomly initialized. Then, after training, each set is used to project the input embeddings (or vectors from lower encoders/decoders) into a different representation subspace.).
   - If we do the same self-attention calculation we outlined above, just eight different times with different weight matrices, we end up with eight different Z matrices. The feed-forward layer is not expecting eight matrices – it’s expecting a single matrix (a vector for each word). So we need a way to condense these eight down into a single matrix.
   - concatenate the Z matrices: We concat the matrices then multiply them by an additional weights matrix WO

One thing that’s missing from the model as we have described it so far is a way to account for the order of the words in the input sequence.
To address this, the transformer adds a vector to each input embedding. These vectors follow a specific pattern that the model learns, which helps it determine the position of each word, or the distance between different words in the sequence. The intuition here is that adding these values to the embeddings provides meaningful distances between the embedding vectors once they’re projected into Q/K/V vectors and during dot-product attention.
The formula for positional encoding is described in the paper (section 3.5).

One detail in the architecture of the encoder that we need to mention before moving on, is that each sub-layer (self-attention, ffnn) in each encoder has a residual connection around it, and is followed by a layer-normalization step.

The encoder start by processing the input sequence. The output of the top encoder is then transformed into a set of attention vectors K and V. These are to be used by each decoder in its “encoder-decoder attention” layer which helps the decoder focus on appropriate places in the input sequence:


The self attention layers in the decoder operate in a slightly different way than the one in the encoder:

In the decoder, the self-attention layer is only allowed to attend to earlier positions in the output sequence. This is done by masking future positions (setting them to -inf) before the softmax step in the self-attention calculation.

The “Encoder-Decoder Attention” layer works just like multiheaded self-attention, except it creates its Queries matrix from the layer below it, and takes the Keys and Values matrix from the output of the encoder stack.

The Final Linear and Softmax Layer

The decoder stack outputs a vector of floats. How do we turn that into a word? That’s the job of the final Linear layer which is followed by a Softmax Layer.
The Linear layer is a simple fully connected neural network that projects the vector produced by the stack of decoders, into a much, much larger vector called a logits vector.
Let’s assume that our model knows 10,000 unique English words (our model’s “output vocabulary”) that it’s learned from its training dataset. This would make the logits vector 10,000 cells wide – each cell corresponding to the score of a unique word. That is how we interpret the output of the model followed by the Linear layer.
The softmax layer then turns those scores into probabilities (all positive, all add up to 1.0). The cell with the highest probability is chosen, and the word associated with it is produced as the output for this time step.

Recap Of Training

Now that we’ve covered the entire forward-pass process through a trained Transformer, it would be useful to glance at the intuition of training the model.
During training, an untrained model would go through the exact same forward pass. But since we are training it on a labeled training dataset, we can compare its output with the actual correct output.

Once we define our output vocabulary, we can use a vector of the same width to indicate each word in our vocabulary. This also known as one-hot encoding
Now, because the model produces the outputs one at a time, we can assume that the model is selecting the word with the highest probability from that probability distribution and throwing away the rest. That’s one way to do it (called greedy decoding). Another way to do it would be to hold on to, say, the top two words (say, ‘I’ and ‘a’ for example), then in the next step, run the model twice: once assuming the first output position was the word ‘I’, and another time assuming the first output position was the word ‘a’, and whichever version produced less error considering both positions #1 and #2 is kept. We repeat this for positions #2 and #3…etc. This method is called “beam search”, where in our example, beam_size was two (meaning that at all times, two partial hypotheses (unfinished translations) are kept in memory), and top_beams is also two (meaning we’ll return two translations). These are both hyperparameters that you can experiment with.


---



































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


