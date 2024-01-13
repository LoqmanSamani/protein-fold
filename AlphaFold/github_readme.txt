You will need a machine running Linux, AlphaFold does not support other operating systems. 
Full installation requires up to 3 TB of disk space to keep genetic databases (SSD storage is recommended) 
and a modern NVIDIA GPU (GPUs with more memory can predict larger protein structures).


# Once the script has finished, you should have the following directory structure:

$DOWNLOAD_DIR/                             # Total: ~ 2.62 TB (download: 556 GB)
    bfd/                                   # ~ 1.8 TB (download: 271.6 GB)
        # 6 files.
    mgnify/                                # ~ 120 GB (download: 67 GB)
        mgy_clusters_2022_05.fa
    params/                                # ~ 5.3 GB (download: 5.3 GB)
        # 5 CASP14 models,
        # 5 pTM models,
        # 5 AlphaFold-Multimer models,
        # LICENSE,
        # = 16 files.
    pdb70/                                 # ~ 56 GB (download: 19.5 GB)
        # 9 files.
    pdb_mmcif/                             # ~ 238 GB (download: 43 GB)
        mmcif_files/
            # About 199,000 .cif files.
        obsolete.dat
    pdb_seqres/                            # ~ 0.2 GB (download: 0.2 GB)
        pdb_seqres.txt
    small_bfd/                             # ~ 17 GB (download: 9.6 GB)
        bfd-first_non_consensus_sequences.fasta
    uniref30/                              # ~ 206 GB (download: 52.5 GB)
        # 7 files.
    uniprot/                               # ~ 105 GB (download: 53 GB)
        uniprot.fasta
    uniref90/                              # ~ 67 GB (download: 34 GB)
        uniref90.fasta



The simplest way to run AlphaFold is using the provided Docker script. 
This was tested on Google Cloud with a machine using the nvidia-gpu-cloud-image with 12 vCPUs, 
85 GB of RAM, a 100 GB boot disk, the databases on an additional 3 TB disk, and an A100 GPU. 

# AlphaFold prediction speed

The table below reports prediction runtimes for proteins of various lengths. 
We only measure unrelaxed structure prediction with three recycles while excluding 
runtimes from MSA and template search. All runtimes are from a single A100 NVIDIA GPU.


No. residues 	Prediction time (s)
100 	              4.9
200 	              7.7
300 	              13
400 	              18
500 	              29
600 	              36
700 	              53
800 	              60
900 	              91
1,000 	              96
1,100 	              140
1,500 	              280
2,000 	              450
2,500 	              969
3,000 	              1,240
3,500 	              2,465
4,000 	              5,660
4,500 	              12,475
5,000 	              18,824


# Examples on how to use AlphaFold in different scenarios

# Folding a monomer

>sequence_name
<SEQUENCE> # fasta format 

python3 docker/run_docker.py \
  --fasta_paths=monomer.fasta \
  --max_template_date=2021-11-01 \
  --model_preset=monomer \
  --data_dir=$DOWNLOAD_DIR \
  --output_dir=/home/user/absolute_path_to_the_output_dir

# Folding a homomer (for example: a homomer with 3 copies of the same sequence <SEQUENCE>)

>sequence_1
<SEQUENCE>
>sequence_2
<SEQUENCE>
>sequence_3
<SEQUENCE>

python3 docker/run_docker.py \
  --fasta_paths=homomer.fasta \
  --max_template_date=2021-11-01 \
  --model_preset=multimer \
  --data_dir=$DOWNLOAD_DIR \
  --output_dir=/home/user/absolute_path_to_the_output_dir

# Folding a heteromer (for example: an A2B3 heteromer, i.e. with 2 copies of <SEQUENCE A> and 3 copies of <SEQUENCE B>.)

>sequence_1
<SEQUENCE A>
>sequence_2
<SEQUENCE A>
>sequence_3
<SEQUENCE B>
>sequence_4
<SEQUENCE B>
>sequence_5
<SEQUENCE B>

python3 docker/run_docker.py \
  --fasta_paths=heteromer.fasta \
  --max_template_date=2021-11-01 \
  --model_preset=multimer \
  --data_dir=$DOWNLOAD_DIR \
  --output_dir=/home/user/absolute_path_to_the_output_dir


# AlphaFold output

The outputs include the computed MSAs, unrelaxed structures, relaxed structures, 
ranked structures, raw model outputs, prediction metadata, and section timings.


# The --output_dir directory will have the following structure:

<target_name>/
    features.pkl  # A pickle file containing the input feature NumPy arrays used by the models to produce the structures.
    ranked_{0,1,2,3,4}.pdb  # A PDB format text file containing the predicted structures, after reordering by model confidence. 
    ranking_debug.json  # A JSON format text file containing the pLDDT values used to perform the model ranking, and a mapping back to the original model names.
    relax_metrics.json  # A JSON format text file containing relax metrics, for instance remaining violations.
    relaxed_model_{1,2,3,4,5}.pdb # A PDB format text file containing the predicted structure, after performing an Amber relaxation procedure on the unrelaxed structure prediction
    result_model_{1,2,3,4,5}.pkl  # A pickle file containing a nested dictionary of the various NumPy arrays directly produced by the model
    timings.json  # A JSON format text file containing the times taken to run each section of the AlphaFold pipeline
    unrelaxed_model_{1,2,3,4,5}.pdb # A PDB format text file containing the predicted structure, exactly as outputted by the model.
    msas/  # A directory containing the files describing the various genetic tool hits that were used to construct the input MSA.
        bfd_uniref_hits.a3m
        mgnify_hits.sto
        uniref90_hits.sto



The provided inference script is optimized for predicting the structure of a single protein, and it will 
compile the neural network to be specialized to exactly the size of the sequence, MSA, and templates. 
For large proteins, the compile time is a negligible fraction of the runtime, 
but it may become more significant for small proteins or if the multi-sequence alignments are already precomputed. 




# If you use the code or data in this package, please cite:

@Article{AlphaFold2021,
  author  = {Jumper, John and Evans, Richard and Pritzel, Alexander and Green, Tim and Figurnov, Michael and Ronneberger, Olaf and Tunyasuvunakool, Kathryn and Bates, Russ and {\v{Z}}{\'\i}dek, Augustin and Potapenko, Anna and Bridgland, Alex and Meyer, Clemens and Kohl, Simon A A and Ballard, Andrew J and Cowie, Andrew and Romera-Paredes, Bernardino and Nikolov, Stanislav and Jain, Rishub and Adler, Jonas and Back, Trevor and Petersen, Stig and Reiman, David and Clancy, Ellen and Zielinski, Michal and Steinegger, Martin and Pacholska, Michalina and Berghammer, Tamas and Bodenstein, Sebastian and Silver, David and Vinyals, Oriol and Senior, Andrew W and Kavukcuoglu, Koray and Kohli, Pushmeet and Hassabis, Demis},
  journal = {Nature},
  title   = {Highly accurate protein structure prediction with {AlphaFold}},
  year    = {2021},
  volume  = {596},
  number  = {7873},
  pages   = {583--589},
  doi     = {10.1038/s41586-021-03819-2}
}


# AlphaFold communicates with and/or references the following separate libraries and packages:

    Abseil
    Biopython
    Chex
    Colab
    Docker
    HH Suite
    HMMER Suite
    Haiku
    Immutabledict
    JAX
    Kalign
    matplotlib
    ML Collections
    NumPy
    OpenMM
    OpenStructure
    pandas
    pymol3d
    SciPy
    Sonnet
    TensorFlow
    Tree
    tqdm




