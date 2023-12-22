## Highly accurate protein structure prediction with AlphaFold

![fold20](https://github.com/LoqmanSamani/protein_sa/blob/systembiology/images/fold20.png)


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


