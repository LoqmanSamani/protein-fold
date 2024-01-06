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

After collecting the data, Input primary sequence and MSA features are embedded to MSA representation (The size of this list is hyperparameter we can set, basically it would be the length of the longest sentence in our training dataset) which is a matrix with the shape of 
(number of sequences * length of the sequence). The structural information are also embedded to form pair representation, a matrix with the shape of (length of the sequence * length of the sequence).

why embedding: embedding means to convert the alphabetic information (A, C, G, T) stores in MSA and also in pair representation to numbers (probabilities)(turning each input word into a vector), which are the input of the neural networks in general.

in training phase 75% of training data comes from genetic search (sequence information) approach and the rest 25% comes from template search (pair representations).

co-evolution: Suppose we have a protein where an amino acid with negative charge (say, glutamate) is near to an amino acid with positive charge (say, lysine), although they are both far away in the amino acid sequence. 
This Coulombic interaction stabilises the structure of the protein. Imagine now that the first amino acid mutates into a positively charged amino acid — in order to preserve this contact, 
the second amino acid will be under evolutionary pressure to mutate into a negatively charged amino acid, otherwise the resulting protein may not be able to fold.



![evoformer](https://github.com/LoqmanSamani/protein_sa/blob/systembiology/%CE%B1_fold/images/evoformer.png)

## 3 . Evoformer (evolutionary transformer?) module

Evoformer is the main step in the alphafold algorithm, it is actually a transformer neural network which takes the MSAs and structural information from the last steps
as input and the output of this network is the prediction of protein structures as a graph inference problem in 3D space in which the edges of the
graph are defined by residues in proximity, it means each node in the graph corresponds to a residue (amino acid), and the edges between nodes represent the spatial relationships between those residues.
The objective of this part is to refine the representations for both the MSA and the pair interactions, but also to iteratively exchange information between them.
The Evoformer network consists of 48 (by default) blocks (iterations), each block has an MSA representation and a pair representation as its input and output and processes them within several layers.

Basically we can devide the evoformer (transformer) into two parts Encoder and Decoder. in general the encoder part used as a block in which the input parameters (in alphafold case inputs are MSA representations)
are processed and the decoder part is a block in which the output parameters are processed(in this case pair representation information are the output).
each encoder or decoder contains two main parts a) self-attention b) feed forward neural network (FFN(x) = max(0, X*W1 + b1 )W2 + b2).
in alphafold in each block there are an encoder and a decoder which share information. The weights in each block of evoformer (neural network parameters are not shared with other blocks it means each block saves the weights for its self just like general transformers.)

a) self-attention 

- axial self-attention in the MSA stack (encoder): this part of evoformer network processes MSA representations in two subparts: 

The first step in calculating self-attention (in both cases row (horizontal) attention with pair bias and column attention) is to create three vectors from each of the encoder’s input vectors (in this case the input are MSA representations)
for each MSA, it creates a Query matrix(q), a Key matrix(k), and a Value matrix(v) and then they will be multiplied by weights matrices in each case (this process called positional encoding which means add a matrix of weights to each input embedding. These matrix follow a specific pattern that the model learns, which helps it determine the position of each word, or the distance between different words in the sequence.
The intuition here is that adding these values to the embeddings give the embedded words in the sequences a positional weight, which gives the words in sequences a positional matter.).
The second step in calculating self-attention is to calculate the attention score (The score determines how much focus to place on other parts of the input sentence as we encode a word at a certain position)
This score is calculated by matrix multiplication of the query matrix with the key matrix.
next the calculated score will divided by the radical 32 (this is a hyperparameter and can be changed, the intuition behind this devision is to lead to having more stable gradients) then the results (scores) will be passed through a softmax operation. Softmax normalizes the scores so they’re all positive and add up to 1.
next each value matrix will be multiplied by the softmax scores to keep intact the values of each residue in the MSA representation.

       a linear transformation to the input column (columns in the MSA representation). This involves a linear layer followed by layer normalization. 
       It produces three sets of projections (q, k, and v), which are associated with query, key and value in attention mechanism.

       Wq * X = Q,  Wk * X = K,  Wv * X = V  , here X is the MSA matrix

       g = sigmoid(Linear(MSA)): Applies a linear transformation to msa followed by the sigmoid activation function. 
       This produces the gating factor g which is used to modulate the attention scores.

       self-attention = softmax((Q * K.T)/ sqrt(32)) * V  
       Computes attention scores using the scaled dot-product attention mechanism. 
       It involves taking the dot product of the query q and key k, scaling it by the square root of the dimensionality 


       ColumnAttention = g . sum(self-attention)
  
       g = sigmoid(LinearNoBias(LinearNorm(z))) applies a linear transformation to the input z(pair representation) with no bias term. This is used for incorporating pair biases in the attention mechanism.(During the per-sequence attention in the MSA, we project additional
       logits from the pair stack to bias the MSA attention. This closes the loop by providing information flow from the pair representation back into the MSA representation, ensuring that the overall Evoformer block is
       able to fully mix information between the pair and MSA representations and prepare for structure generation within the structure module.)

       self-attention = softmax((Q * K.T)/ sqrt(32) + b) * V

       RowAttentionWithPairWise = g . sum(self-attention)
       (The most important feature of AlphaFold 2’s MSA transformer is that the row-wise (horizontal) attention mechanism incorporates information from the “pair representation”. When computing attention, the network adds a bias term that is calculated directly from the current pair representation. This trick augments the attention mechanism and allows it to pinpoint interacting pairs of residues.)

![evoformer1](https://github.com/LoqmanSamani/protein_sa/blob/systembiology/%CE%B1_fold/images/evoformer1.png)

![evoformer2](https://github.com/LoqmanSamani/protein_sa/blob/systembiology/%CE%B1_fold/images/evoformer2.png)


After row-wise and column-wise attention the MSA stack contains a 2-layer MLP(multi-layer perceptron) as the transition layer.This stage of processing in the evoformer block operates across features,
refining the representation using a non-linear transform.

![evoformer3](https://github.com/LoqmanSamani/protein_sa/blob/systembiology/%CE%B1_fold/images/evoformer3.png)

The “Outer product mean” block transforms the MSA representation into an update for the pair representation . All MSA entries are linearly projected to a smaller dimension c = 32 with
two independent Linear transforms. The outer products of these vectors (If vi and vj are the vectors obtained from the linear projections for columns i and j, then the outer product is vi⊗vj)from two columns i and j are averaged over the sequences (This involves taking the mean over the corresponding elements of the outer products across the sequences.) and projected to dimension cz to obtain an update for entry ij in the pair representation.
Mathematically, if viand vj are column vectors obtained from linear projections, and ⊗ denotes the outer product, the update for entry ij (Uij) can be expressed as:

      Uij = Pij(mean(vi⊗vj))
Here, Pij is a linear transform, and mean calculates the mean over the sequences.

![evofomer4](https://github.com/LoqmanSamani/protein_sa/blob/systembiology/%CE%B1_fold/images/evoformer4.png)

next step in evoformer block is The triangular multiplicative update (The key feature of this network is that attention is arranged in terms of triangles of residues. 
The intuition here is to enforce the triangle inequality, one of the axioms of metric spaces. This is quite a clever idea since one of the classical problems of deep learning-based 
structure prediction was that distance distributions could not be embedded in three-dimensional space) updates the pair representation in the Evoformer block by
combining information within each triangle of graph edges ij, ik, and jk. Each edge ij receives an update
from the other two edges of all triangles, where it is involved. this step sontains two sub steps, one for the
“outgoing” edges and one for the “incoming” edges.
- The input pair representation Zij is first normalized
- Two linear projections are applied to Zij to obtain Aij and Bij

     Aij,Bij = sigmoid(Linear(Zij))⋅Linear(Zij)
     Linear denotes the linear transformation, and sigmoid is the sigmoid activation function.

- Another linear transformation Gij is applied to Zij and passed through a sigmoid activation: 

     Gij = sigmoid(Linear(Zij))

- The multiplicative update is performed using the computed Gij and a linear transformation of the sum of element-wise products:

    Z_hat ij=Gij⋅Linear(LayerNorm(∑(Aik ⋅Bjk)))
    The result is the updated pair representation Z_hat ij

![evoformer5&6](https://github.com/LoqmanSamani/protein_sa/blob/systembiology/%CE%B1_fold/images/evoformer5%266.png)


next step in evoformer block is Triangular self-attention: this process takes place in two sub steps:
first the “starting node” version updates the edge ij with values from all edges that share the same starting node i, (i.e. all edges ik). The decision whether edge ij will receive an update from edge ik is
not only determined by their query-key similarity, but also modulated by the
information bjk derived from the third edge jk of this triangle. the second sub step updates with an
additional gating gij derived from edge ij. The symmetric pair of this module operates on the edges around
the ending node.
here the same steps occur as in the self-attention in MSA representations.

![evoformer7]()


next step is The transition layer. this step in the pair stack is equivalent to that in the MSA stack: a 2-layer MLP
where the intermediate number of channels expands the original number of channels by a factor of 4.


#### importaint information about evoformer:

- The objective of attention is to identify which parts of the input are more important for the objective of the neural network. In other words, to identify which parts of the input it should pay attention to.
- The final Evoformer block provides a highly processed MSA representation {msi } and a pair representation {zij }, which contain information required for the structure module.





## 4. The structure module

![structure](https://github.com/LoqmanSamani/protein_sa/blob/systembiology/%CE%B1_fold/images/structure.png)

this part of the alphafold algorithm takes the refined MSA representation and pair representation, which were manipulated with Evoformer, and leverages them to construct a three-dimensional model of the structure.
The end result is a long list of Cartesian coordinates representing the position of each atom of the protein, including side chains.

the evoformer single representation (The prediction modules are also using a “single” sequence representation {si} with si ∈ Rcs , cs = 384 and i ∈ {1 . . . Nres }.)
is used as the initial single representation and Evoformer’s pair representation ({zij } with zij ∈ Rcz and i, j ∈ {1, ..., Nres }) biases the affinity maps in the attention operations.
The module has 8 layers with shared weights. Each layer updates the abstract single representation {si } as
well the concrete 3D representation (the “residue gas”) which is encoded as one backbone frame per residue
{Ti=(Ri, ti)}.
The Ti represents an Euclidean transform from the local frame to a global reference frame. I.e. it transforms a position in local coordinates ~X_local ∈ R3 to a position in global coordinates ~X_global ∈ R3
      
     ~X_global = Ti ◦ ~X_local
              = Ri~X_local + ~ti

The backbone frames are initialized to an identity transform. We call this approach
the ‘black hole initialization’. This initialization means that at the start of the Structure module all residues
are located at the same point (the origin of the global frame) with the same orientation.

One “layer” of the structure module is composed by the following operations:




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


