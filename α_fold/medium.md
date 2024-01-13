## Highly accurate protein structure prediction with AlphaFold

### Introduction
Proteins, which are fundamental to life, play a crucial role in biological processes, so understanding their structures is essential for deciphering their functions. 
Despite experimental efforts that have led to structural insights for about 100,000 unique proteins, this is only a fraction of the millions of known protein sequences. 
The bottleneck in structural coverage is due to the laborious and time-consuming nature of experimental determination. 
To close this gap and enable large-scale structural bioinformatics, accurate computational approaches are essential. Predicting the three-dimensional structure of a protein 
based solely on its amino acid sequence, a key aspect of the protein folding problem, has been a challenge for over 50 years. 
While experimental structures cover only 17% of human protein residues, AlphaFold2 represents a breakthrough computational solution that utilizes machine learning methods 
on an unprecedented scale to predict protein structures with atomic accuracy. 

AlphaFold covers a remarkable 98.5% of the human proteome and provides reliable predictions for 58% of residues, with an exceptional subset showing a very high degree of confidence.
In this article, I would like to explain the intricate details of the AlphaFold system (Figure 1) and break down its methodology step by step. As illustrated in Figure 1 the system can be divided into 
three main components: **Database Search**, **Evoformer Module** and **Structure Module**. 
Each of these segments plays a pivotal role in the system, and the flow of information retrieved from the protein database through the various parts of the system leads to an extremely 
accurate 3D structure of a protein.


![alphafold](https://github.com/LoqmanSamani/protein_sa/blob/systembiology/%CE%B1_fold/images/alphafold.png)
<figcaption style="font-size: smaller; text-align: center;">Figure 1: The AlphaFold system. The system consists of three main components: Database Search, Evoformer Module and Structure Module. The detailed explanations of each component and its layers (represented by numbers) are detailed in the corresponding part of the article.</figcaption>



### Genetic and Structure Database Searches

The initial stage of the AlphaFold system is the **Data Pipeline**, which is the process of gathering both sequential and structural information essential for effective model training and enhancing the
accuracy of the inference process.
In addition to amino acid sequence files in the **FASTA format**, another resource utilized during the training phase of the system is the **Macromolecular Crystallographic Information File (mmCIF)**. 
As implied by its name, mmCIF contains a diversity of sequential and structural protein information derived from proteins collected through experimental methods such as **X-ray crystallography**. 
During the inference phase, the system requires only the sequence identifier, which is then used to retrieve the amino acid sequence of the protein in the FASTA format. 
The system parses the input files (mmCIF and fasta) and extracts relevant information, in the case of the mmCIF this is the **name of the sequence**, the **amino acid sequence**, 
the **atomic coordinates**, the **release date** and the **resolution**. In the other case (fasta), the amino acid sequence and its name are the only information required.
With the information retrieved from the files, two different database searches are performed, in one case, the **genetic database search** (Figure 1. 2), the system searches
**the Big Fantastic Database (BFD)** for similar sequences (in this case, the system used specialized applications such as **JackHMMER** and **HHBlits**) and the result is called 
**multiple sequence alignment(MSA)**, then the result is modified by removing duplicate sequences and those that are not long enough to be included in the process.
The remaining sequences are then assembled into a matrix of sequences called the *MSA representation*. Each column of this matrix encodes a single residue of the 
input sequence in different organisms, while each row represents the entire sequence in a particular organism.
In another case, the system uses an approach called **HHSearch** to search in the **Protein Data Bank (PDB)**, **Structure Database Search** (Figure 1. 3), to retrieve the 3d structures of the proteins, 
that have a similar primary structure to the original input(Figure 1. 1). 

### Embedding process

Before the collected sequential (MSA information) and structural information can be used in the alphafold system, it must be embedded. 
This is the process of converting categorical information, such as the alphabetical representations (G, V, A, ... . Each letter represents an amino acid) 
found in MSA, as well as the structural features, into numerical vectors or probabilities. 
For this purpose, each input symbol or word is converted into a vectorized representation. These numerical vectors serve as input for neural networks, 
enabling them to effectively process and learn patterns from the sequential data.
Regarding the MSA, the resulting embedded information is denoted as the **MSA representation** (Figure 1. 4). 
This representation takes the form of a matrix with dimensions (number of sequences * length of the longest sequence). 
Conversely, in the other case, the embedded outcome is termed the **Pair Representation** (Figure 1. 5), which manifests as a matrix with dimensions (length of sequence * length of sequence). 
These distinct yet crucial matrices serve as the primary input parameters during both the training and inference phases of the system. 
They encapsulate essential structural and sequential information, facilitating subsequent steps in the algorithm.


### Evoformer module

The Evoformer (Figure 1. 6) stands out as a pivotal and potent module within the AlphaFold system. Fundamentally, it operates as a variant of the transformer neural network, originally introduced by Google Brain in 2017. 
Initially designed for language translation tasks, transformers excel with sequential data. In the context of AlphaFold, the Evoformer takes both MSA-representation and pair representation as input. 
Leveraging this information, it produces an output that is a prediction of the protein structure as a graph inference problem in 3D space. In this 3D space, 
the nodes of the graph correspond to amino acid residues in proximity, while the edges represent the spatial relationships between these residues. 

![evoformer](https://github.com/LoqmanSamani/protein_sa/blob/systembiology/%CE%B1_fold/images/evoformer.png)
<figcaption style="font-size: smaller; text-align: center;">Figure 2: The Evoformer network. it consists of 48 (by default) blocks, each block has an MSA representation(1) and a pair representation(2) as its input and output and processes them within several layers. The detailed explanations of each layer, the core components of this network's functionality, are meticulously expounded in the corresponding part of the article.</figcaption>

The Evoformer block (figure 2) comprising two transformer models, each dedicated to processing one of the primary input datasets, MSA representation (Figure 2, 1) and pair representation (Figure 2, 2). 
These transformers collaboratively exchange information across two distinct layers within the block. 
In the initial layer of the block, information flows from pair representation to MSA representation. Meanwhile, in the fourth layer (outer product mean), information is reciprocally 
transmitted from the MSA representation transformer to the pair representation transformer. This intricate information exchange mechanism enhances the Evoformer's ability to 
integrate insights from both MSA and pair representations, contributing to the overall refinement and accuracy of the protein structure prediction.
Within each of the two transformers in the Evoformer block, there are two main components. The initial part is referred to as **self-attention**, 
which allows the model to weigh the significance of different positions within the input sequence. The second part is denoted as the **Feed Forward Neural Network (FFNN)**, 
which processes the information from the self-attention mechanism, aiding in capturing intricate patterns and dependencies within the data. 

#### Self Attention

The main idea behind self-attention in a Transformer is to enable the model to weigh the importance of different words in a sequence when processing each word. 
It allows the model to consider the relationships and dependencies between words in a flexible and adaptive way, capturing long-range dependencies and improving 
its ability to understand the context of each word.
In the context of protein structure, a diverse set of 20 unique amino acids (analogous to words) in each sequence plays a vital role in determining distinct folding patterns. 
It is essential to gather detailed information about each residue and understand how they interact with neighboring residues. 
This comprehensive knowledge is crucial for inferring higher-order structures such as secondary, tertiary, and quaternary structures in proteins.
To gain this necessary information alphafold defines **axial self-attention** (Figure 2. 3 & 4)in the MSA representation transformer and **triangular self-attention** (Figure 2. 9 & 10) in the pair representation transformer.

##### axial self-attention in the MSA stack

The initial step in the Evoformer's self-attention mechanism involves creating three vectors—query (q), key (k), and value (v)—from each row (row-attention, Figure 2. 3) or column (column-attention, Figure 2. 4) 
of the embedded MSA-representation. For each vector, a corresponding weight matrix (Wq, Wk, and Wv) is employed, and these matrices are then multiplied by the embedded MSA (X) to generate 
the Q, K, and V matrices. This procedure is known as positional encoding, wherein a matrix of weights is added to each input embedding. The positional encoding matrix adheres to a specific 
pattern learned by the model, aiding in determining the position of each residue or the distance between different residues in the sequence.

Following the creation of the Q, K, and V matrices, the next step in the Evoformer's self-attention mechanism involves normalizing these matrices using the softmax function. 
This normalization ensures that the scores are positive and sum up to 1. As depicted in Figure 2-a, there exists a connection between the MSA-transformer and pair-transformer,
which demonstrated the flow of information from pair representation to the MSA-representation. In this process, additional information from pair representation is incorporated
to influence the bias terms of the MSA attention.
    
    Wq * X = Q
    Wk * X = K
    Wv * X = V

    column-attention: Z = softmax((Q . K.T)/sqrt(d)) * V;     d = dimension of the keys, queries, and values matrices
    row-attention: Z = softmax((Q . K.T)/sqrt(d) + b) * V;    b = pair representation information

#### Feed Forward Neural Networks

After row-wise and column-wise attention next layer of the MSA transformer contains a 2-layer **MLP, multi-layer perceptron**, known as FFNN (Figure 2. 5),  as the transition layer.This stage of processing in the Evoformer block operates across features,
refining the representation using a non-linear transform. The main idea behind a feed-forward neural network is to process input data through a series of layers, where each layer consists of nodes, connected to nodes in the subsequent layer. 
During the training process of a Feedforward Neural Network (FFNN), each weight matrix associated with the layers of the network is optimized using **Adam optimization algorithm** to improve the accuracy of predicted structures. Notably, 
in a conventional neural network architecture, weights are typically shared across different layers. 
However, in the case of the Evoformer, each block operates with its own set of weights, and these weights are not shared with other blocks. 
This contrasts with the usual neural network setup, where weights are often shared across layers. The unique characteristic of Evoformer lies in the independence of weights for each block, 
allowing for more localized and specific adaptations during the training process.


    FFN(X) = max(0, X * W1 + b1 ) * W2 + b2


#### Outer Product Mean

This layer of the Evoformer (Figure 2. 6) serves the purpose of transforming the MSA-representation into an update for the pair representation. 
The transformation involves computing the mean of the outer product of each two columns (for instance, Ci and Cj), and this result is utilized as an update for the element (Cij) in the pair representation. 
This mechanism ensures the integration of information from the MSA into the Pair representation module.

Mathematically, this process can be represented as:

     outer product matrix = C¹(m*n) * C²(m*n).T  
     outer product mean = C³ = mean(outer product mean)


#### Triangle multiplicative update and triangle self-attention

In the Evoformer block, the process of updating the pair representation is designed to implement an attention mechanism akin to the one applied in the Multiple 
Sequence Alignment (MSA) representation. The objective is to capture intricate relationships between residues in a three-dimensional protein structure, enabling 
a nuanced understanding of how each residue influences the spatial representation of others in the molecular space.

As illustrated in Figure 3a, a graph representation is constructed from the pair representation matrix. In this graph, each node represents an amino acid, 
and the edges correspond to entries in the matrix. The circles within the graph symbolize individual residues. To simplify the explanation, we focus on 
updating a single edge, ij, leveraging information from other edges.

During the Triangular Multiplicative Update step, the edge ij undergoes an update by assimilating information from the other edges. 
It is noteworthy that there exist two symmetric versions of this update, differentiating between **outgoing** and **incoming** edges. 
This ensures a comprehensive consideration of information flow within the graph structure.

Following the Triangular Multiplicative Update, the process continues with the Triangular Self-Attention step. 
In this step, the edge ij is further updated by incorporating values from all edges that share the same starting node. 
This comprehensive attention mechanism facilitates the integration of information from neighboring residues, 
contributing to a refined and context-aware representation of the 3D spatial relationships between amino acids in the protein structure.

After the attention mechanism processes the pair representation in the Evoformer block, the subsequent step involves a transition layer (feedforward neural network). 
This transition layer is analogous to the one found in the MSA representation transformer. In this layer, the pair representation is further refined and transformed 
through the application of a feedforward neural network.

![structure_attention](https://github.com/LoqmanSamani/protein_sa/blob/systembiology/%CE%B1_fold/images/structure_attention.png)
<figcaption style="font-size: smaller; text-align: center;">Figure 3: Triangle Multiplicative Update and Triangle Self-Attention. a) The pair representation interpreted as directed edges in a graph. b) Triangle multiplicative update and triangle self-attention. The circles represent residues. Entries in the pair representation are illustrated as directed edges, and in each diagram, the edge being updated is ij.</figcaption>




### The structure module

![structure](https://github.com/LoqmanSamani/protein_sa/blob/systembiology/%CE%B1_fold/images/structure.png)

In the final step of AlphaFold system, the algorithm extracts the first row of the highly processed MSA representation, which corresponds to the original sequence, and utilizes the pair representation. 
These two representations serve as the foundation for the construction of a three-dimensional model of the protein structure (structure 6)
The outcome of this process is a comprehensive list of **Cartesian coordinates** presented in the Protein Data Bank (PDB) format. These coordinates delineate the precise positions of 
each atom within the protein structure, encompassing not only the backbone but also the intricate details of side chains. 



the single representation (structure 1) and pair representation and backbone frames(spatial representations)
are used as input of the structure module. This module works with a gradient descent base algorithm to optimize the predicted 3d structure (spatial representations).
This module utilizes a gradient descent-based algorithm to optimize the predicted 3D structure of the protein. The structure module comprises 8 iterations by default, 
and in each iteration, the weights are shared among the layers, differing from the Evoformer module where each block has distinct weights.

During each iteration of the Structure Module, the following steps take place:

### Invariant Point Attention
In this step of the Structure Module, known as **Invariant Point Attention (IPA)**, a specialized type of attention mechanism designed for handling 3D structures is employed. 
The general attention mechanism, previously explained in detail in the Evoformer module, is adapted for use in this specific context.
The Pair Representation and Single Representation (first row of the MSA representation) play important roles in this sub-step. These representations are utilized to iteratively refine the backbone frames, 
which collectively define the 3D structure of the protein.

In each iteration of the Structure Module during the training phase, the information (backbone frames) generated by IPA and the Single Representation is utilized to predict the 3D structure of the protein. 
The predicted structure at the end of each iteration is then compared to the actual structure using two key metrics:

#### Frame Aligned Point Score

This metric measures the difference between all atom predicted coordinates and the actual coordinates.
It provides  assessment of the alignment between the predicted and actual positions of each atom in the protein structure.

#### Torsion Angle Loss:
        
**Torsion angles** represent the rotations around the chemical bonds in the protein structure, 
including both side chain and backbone torsion angles. The torsion angle loss quantifies the difference between the predicted torsion angles and the actual angles.

The calculated values for Frame Aligned Point Score and Torsion Angle Loss are combined, to form a unified **loss function**. This combined loss is then employed in the gradient descent algorithm during the training process.
The gradient descent algorithm adjusts the model's parameters to minimize this loss, facilitating the refinement of the predicted 3D protein structure. 
By iteratively optimizing the model based on the disparity between predictions and actual observations, the algorithm converges towards a more accurate and biologically relevant protein structure.



### Recycling

Following the prediction of the final protein structure in the Structure Module, all the generated and refined information, including the Multiple Sequence Alignment (MSA) representation, 
pair representation, and the predicted 3D structure, serves as input for Evoformer blocks. This recycling mechanism enables the model to iteratively refine its predictions.
This iterative refinement process is repeated three times by default. In each iteration, the output from the previous step is embedded and utilized as additional input to the model. 
This iterative embedding of information ensures that the model progressively refines its understanding and representation of the protein structure. 
By incorporating the refined predictions from earlier iterations, the model can capture finer details, improve precision, 
and enhance the overall accuracy of its final predictions. The recycling mechanism contributes to the robustness and effectiveness 
of the AlphaFold algorithm in predicting highly accurate protein structures.



### Conclusion

In conclusion, AlphaFold represents a groundbreaking advance in the field of protein structure prediction, bridging the gap between experimental efforts and the vast landscape of known protein sequences. 
The intricate interplay of its three main components - database search, Evoformer module and structure module - contributes to the remarkable achievement of predicting protein structures with unprecedented 
accuracy. Using machine learning methods to an impressive extent, AlphaFold not only covers 98.5% of the human proteome, but also provides reliable predictions for 58% of residues, 
with a subset showing an exceptional level of confidence. This groundbreaking computational approach has the potential to significantly impact structural bioinformatics by overcoming 
the challenges posed by laborious experimental determinations and providing insights into the three-dimensional shape of proteins that are critical to understanding their functions. 
As we delve into the details of the AlphaFold methodology, it becomes clear that the fusion of genetic database searches, Evoformer's neural network architecture and the structural module 
results in an exceptionally accurate 3D representation of proteins. The success of AlphaFold underscores the power of computational methods in unlocking the secrets of protein structures, 
providing insight into the intricate world of molecular biology and opening new avenues for research and discovery.





### References

[1] Jumper, J., Evans, R., Pritzel, A. et al. [Highly accurate protein structure prediction with AlphaFold](https://www.nature.com/articles/s41586-021-03819-2#citeas). Nature 596, 583–589 (2021)

[2] Tunyasuvunakool, K., Adler, J., Wu, Z. et al. [Highly accurate protein structure prediction for the human proteome](https://www.nature.com/articles/s41586-021-03828-1#citeas). Nature 596, 590–596 (2021).

[3] Senior, A.W., Evans, R., Jumper, J. et al. [Improved protein structure prediction using potentials from deep learning](https://www.nature.com/articles/s41586-019-1923-7#citeas). Nature 577, 706–710 (2020)

[4] Kiersten M. Ruff and Rohit V. Pappu [AlphaFold and Implications for Intrinsically Disordered Proteins](https://www.sciencedirect.com/science/article/pii/S0022283621004411). Journal of Molecular Biology, Volume 433, Issue 20, (2021).

[5] Pearce R, Zhang Y. [Deep learning techniques have significantly impacted protein structure prediction and protein design](https://pubmed.ncbi.nlm.nih.gov/33639355/). Curr Opin Struct Biol. 2021 Jun;68:194-207. doi: 10.1016/j.sbi.2021.01.007. Epub 2021 Feb 24. PMID: 33639355; PMCID: PMC8222070.

[6] Senior AW, Evans R, Jumper J, Kirkpatrick J, Sifre L, Green T, Qin C, Žídek A, Nelson AWR, Bridgland A, Penedones H, Petersen S, Simonyan K, Crossan S, Kohli P, Jones DT, Silver D, Kavukcuoglu K, Hassabis D. [Protein structure prediction using multiple deep neural networks in the 13th Critical Assessment of Protein Structure Prediction (CASP13)](https://pubmed.ncbi.nlm.nih.gov/31602685/). Proteins. 2019 Dec;87(12):1141-1148. doi: 10.1002/prot.25834. PMID: 31602685; PMCID: PMC7079254.

[7] Ashish Vaswani, Noam Shazeer, Niki Parmar, Jakob Uszkoreit, Llion Jones, Aidan N. Gomez, Lukasz Kaiser, Illia Polosukhin. [Attention Is All You Need](https://arxiv.org/abs/1706.03762)

[8] Richard E. Turner. [An Introduction to Transformers](https://arxiv.org/abs/2304.10557)







------------------------------------------------------
------------------------------------------------------
------------------------------------------------------
------------------------------------------------------
------------------------------------------------------






































mathematically equivalent to the cosine of the angle difference (see below).
![ipa](https://github.com/LoqmanSamani/protein_sa/blob/systembiology/%CE%B1_fold/images/ipa.png)

   first for each token(residue) in single representation  the query, key, and value vectors will be created (the same as the general attention mechanism)
   to contain information from the 3d representation and have a mechanism by which the single representation will be affected by spatial information the IPA used the following equation:
![ipa1](https://github.com/LoqmanSamani/protein_sa/blob/systembiology/%CE%B1_fold/images/ipa1.png)   
 
![ipa4](https://github.com/LoqmanSamani/protein_sa/blob/systembiology/%CE%B1_fold/images/ipa4.png)

2) The Structure Module predicts backbone frames Ti and torsion angles if . The atom coordinates are then constructed by applying the torsion angles to the corresponding amino acid structure with idealized bond angles and bond lengths.
in each iteration of structure module In order to resolve any remaining structural violations and clashes, we relax our model predictions by
an iterative restrained energy minimization procedure. At each round, we perform minimization of the AMBER99SB[100] force field with additional harmonic restraints that keep the system near its input structure.

at the end of the structural module this accuracy numbers will be calculated:
1) Frame aligned point error (FAPE): scores a set of predicted atom coordinates under a set of predicted local frames against the corresponding ground truth atom coordinates and ground truth local frames.
   The final FAPE loss scores all atoms in all backbone and side chain frames.
2) Model confidence prediction (pLDDT): compute the expected value of the per-residue pLDDT distribution
3) TM-score prediction: The pLDDT head above predicts the value of lDDT-Cα, which is a local error metric that operates pairwise
   but by design is not sensitive to what fraction of the residues can be aligned using a single global rotation
   and translation. This can be disadvantageous for assessing whether the model is confident in its overall
   domain packing for large chains. In this section we develop a predictor of the global superposition metric
   TM-score.

   
![table1](https://github.com/LoqmanSamani/protein_sa/blob/systembiology/%CE%B1_fold/images/table.png)








After updating pair representation by MSA representation in the next three steps pair representation will be processed.
The triangular multiplicative update updates the pair representation in the Evoformer block by
combining information within each triangle of graph edges ij, ik, and jk. Each edge ij receives an update
from the other two edges of all triangles, where it is involved. There are two symmetric versions, one for the
“outgoing” edges and one for the “incoming” edges.

next step in evoformer block is Triangular self-attention: this process takes place in two sub steps:
first the “starting node” version updates the edge ij with values from all edges that share the same starting node i, (i.e. all edges ik). 
 the second sub step updates with an
additional gating gij derived from edge ij. The symmetric pair of this module operates on the edges around
the ending node.
here the same steps occur as in the self-attention in MSA representations.






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




     







The result of self-attention step (matrix Z) captures attention information which will be sent to the feed forward neural network to

in both cases row-attention (figure 2, a) and column-attention (figure 2, b)

each encoder or decoder contains two main parts a) self-attention b) feed forward neural network (FFN(x) = max(0, X*W1 + b1 )W2 + b2
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

       self-attention = softmax((Q * K.T)/ sqrt(d) + b) * V   # d = dimension of the keys, queries, and values matrices

       RowAttentionWithPairWise = g . sum(self-attention)
       (The most important feature of AlphaFold 2’s MSA transformer is that the row-wise (horizontal) attention mechanism incorporates information from the “pair representation”. When computing attention, the network adds a bias term that is calculated directly from the current pair representation. This trick augments the attention mechanism and allows it to pinpoint interacting pairs of residues.)

![evoformer1](https://github.com/LoqmanSamani/protein_sa/blob/systembiology/%CE%B1_fold/images/evoformer1.png)

![evoformer2](https://github.com/LoqmanSamani/protein_sa/blob/systembiology/%CE%B1_fold/images/evoformer2.png)




![evoformer3](https://github.com/LoqmanSamani/protein_sa/blob/systembiology/%CE%B1_fold/images/evoformer3.png)


![evofomer4](https://github.com/LoqmanSamani/protein_sa/blob/systembiology/%CE%B1_fold/images/evoformer4.png)


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









--------------------------------------------------------------------------------------------------------------------------------------------------



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

Exact enforcement of peptide bond geometry is only achieved in the post-prediction relaxation of the structure by gradient descent in the Amber32 force field.

The residue gas representation is updated iteratively in two stages. 1) a geometry-aware attention operation that we term ‘invariant point attention’ (IPA) is used to update an Nres set of neural
activations (single representation) without changing the 3D positions, 2) an equivariant update operation is performed on the residue gas using the updated activations.

The IPA augments each of the usual attention queries, keys and values with 3D points that are produced in the local frame of each residue such that the final value is invariant to global rotations and translations

The IPA module combines the pair representation, the single representation and the geometric representation to update the single representation
The IPA operates in 3D space. Each residue produces query points, key points and value points in its local frame. These points are projected into the global frame using the backbone frame of the residue in which they interact with each
other. The resulting points are then projected back into the local frame. The affinity computation in the 3D space uses squared distances and the coordinate transformations ensure the invariance of this module
with respect to the global frame.


----



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


co-evolution: Suppose we have a protein where an amino acid with negative charge (say, glutamate) is near to an amino acid with positive charge (say, lysine), although they are both far away in the amino acid sequence. 
This Coulombic interaction stabilises the structure of the protein. Imagine now that the first amino acid mutates into a positively charged amino acid — in order to preserve this contact, 
the second amino acid will be under evolutionary pressure to mutate into a negatively charged amino acid, otherwise the resulting protein may not be able to fold.

