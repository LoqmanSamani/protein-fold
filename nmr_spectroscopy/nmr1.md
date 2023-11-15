##  Nuclear Magnetic Resonance (NMR) Spectroscopy In Biology


**The process of resonance assignment:** This step of NMR is compared with a jigsaw puzzle

If the 3D-Structure of protein is already known, NMR can be used to study the interactions between molecules, providing information about binding interfaces, kinetics, and thermodynamics.

NMR can also be employed to monitor the dynamic movements of proteins, which are critical for their biological function.

### NMR Evolution

-  NMR was independently discovered by two groups led by Felix Bloch and Edward Purcell at the end of World War II.
   (the Nobel Prize in Physics in 1952).

-  This first observation of the **chemical shift**: Proctor and Yu.

-  1952, producing the first commercial **NMR spectrometer**.

-  1953, Overhauser observed the **nuclear Overhauser effect**.

-  1960s, introduction of **superconducting magnets**.

-  1966, **Fourier transformed (FT)** NMR.

-   1976, introduction of **two-dimensional (2D) NMR** by R. Ernst and others. 

-  1985, the first structure of a small globular protein.

-  1991, Nobel Prize in Chemistry for  Richard Ernst. **high-resolution NMR spectroscopy**.

-  2002, Nobel Prize in Chemistry for Kurt Wüthrich. The three-dimensional structures of biological macromolecules in solution.

- 

### 2D and multi-dimensional NMR

In the early days of NMR, 2D NMR was introduced for studying biomolecules by K. Wüthrich's group. They demonstrated its feasibility using a 10 mM sample of BPTI, a 58 amino acid protein, despite limited computational resources at the time. However, high protein concentrations were initially required, which limited its practical use. Over the years, technical improvements reduced this constraint.

***One-Dimensional NMR (1D NMR)*** : A standard NMR spectrum is typically a one-dimensional plot, showing the absorption of nuclear spins at specific resonance frequencies. In 1D NMR, each nucleus has its own signal at a particular frequency, typically represented in parts per million (ppm). This method offers limited information about individual nuclei.

#### 2D NMR:

In the early 1970s, 2D NMR was introduced as a revolutionary technique to correlate the resonance frequencies of different nuclei. Unlike 1D NMR, 2D NMR experiments provide additional information, including the nature of the interacting partners (such as J-coupling or nuclear Overhauser effects) and the interaction strengths. 

***Types of 2D NMR Experiments :***     ***COSY (Correlation Spectroscopy)*** Correlates 1H-1H J-couplings and reveals connectivities between protons in a molecule.
***TOCSY (Total Correlation Spectroscopy)*** Provides a detailed view of spin systems and can distinguish between coupled and uncoupled protons.
***NOESY (Nuclear Overhauser Effect Spectroscopy)*** Reveals nuclear Overhauser effects, which are indicative of spatial proximity between nuclei.
***HSQC (Heteronuclear Single Quantum Correlation)*** Correlates 1H with another nucleus like 15N or 13C, offering valuable structural information.

***Application of 2D NMR*** : 2D NMR spectra, such as 1H-15N HSQC, are often the first NMR experiments recorded for protein structural analysis. These spectra can confirm protein compactness and provide information about the aggregation state of the protein. 2D NMR experiments are robust and require only small amounts of labeled material.

***Multi-Dimensional NMR*** : The modular design of 2D NMR can be extended to 3D or even 4D NMR experiments. Multi-dimensional NMR experiments correlate more nuclei and provide higher-resolution data. These experiments are particularly valuable for structural and dynamic studies of complex biomolecules like proteins.

### Chemical Shift

It refers to the position of a resonance signal in an NMR spectrum and is typically expressed in parts per million (ppm) relative to a reference compound's signal frequency (tetramethylsilane (TMS) is a common reference in organic solvents. In protein NMR, 2,2-dimethyl-2-silapentane-5-sulfonic acid (DSS) is preferred. Some researchers use the water signal as a reference, although its position can be influenced by temperature and pH). This allows for easy comparison of chemical shifts between NMR data recorded at different magnetic field strengths.

***Interpreting Chemical Shifts*** :  If the chemical shifts of one compound change when another is added to the sample, it indicates an interaction. When the resonances of the compound are assigned, these changes can be interpreted at the atomic level. This is particularly useful in understanding protein-ligand interactions, helping identify the interacting parts of the small molecule and the target protein.

***Anisotropic Interaction***: Chemical shift is inherently an anisotropic interaction, but in solution NMR, only the isotropic part is observed. High magnetic fields can induce chemical shift anisotropy (CSA), causing broadening of NMR signals for some nuclei.

***Sensitivity to Structural Effects***: Chemical shifts are highly sensitive to steric and electronic effects and, in the case of proteins, secondary and tertiary structure. Unlike nuclear Overhauser effects (nOe) and J-coupling, chemical shift doesn't depend on specific pairwise interactions between well-identified partners, making its prediction and interpretation more complex.

***Databases and Insights***: As more proteins are assigned in NMR studies, researchers have gained insights into the factors affecting chemical shifts, such as torsion angles, aromatic rings, solvent exposure, temperature, pH, and ionic strength. There are several online databases, with the BioMag-ResBank being the largest, containing thousands of entries. Smaller curated databases like TALOS or TALOS+ have also been developed for specific research purposes. These databases provide valuable reference data for interpreting chemical shifts in NMR spectra.

The ***Chemical Shift Index*** (CSI) is a method that uses chemical shifts to identify the type and location of secondary structures along a protein chain. It can provide information at the residue level, allowing for the determination of secondary structure without the need for additional measurements or structural computations

***Torsion Angles Prediction***: Chemical shifts can also be used to predict torsion angles, which define the backbone conformation of amino acids. 

***Protein Structure Prediction*** : Ongoing research is focused on computing protein structures using only chemical shift information, which can simplify the process and save time compared to traditional methods that rely on nuclear Overhauser effect (nOe) measurements. The CS-Rosetta approach uses Rosetta, an algorithm for de novo protein modeling, and incorporates backbone chemical shifts to select suitable fragments. This approach is particularly valuable for proteins with relatively simple topologies and has the potential to streamline protein structure determination.


### Scalar Coupling (J-Coupling): 

Scalar coupling is an indirect interaction between two nuclei (labeled as A and X) with non-zero spins in a molecule. This interaction occurs through the shared electrons in the molecule. One nucleus perturbs the electron spins, which, in turn, perturb the other nucleus. 

Scalar coupling causes NMR signals to be split into multiple peaks. For example, if two spins (A and X) are scalar coupled, the NMR spectrum of each of these spins will appear as a doublet, as shown in Fig. 4. The separation between these two lines is known as the coupling constant (JAX). This splitting occurs because the spins of A are influenced by the spin state of their neighbor X.

### The Nuclear Overhauser Effect (nOe) 

The dipolar coupling can be observed in solid-state NMR or when the molecular environment is no longer isotropic, achieved by introducing specific alignment media like lipid bicelles or bacteriophages.

The nOe is a short-range through-space interaction observed in NMR. It refers to intensity changes in the NMR spectrum of one nucleus (e.g., spin 1) when another nucleus (e.g., spin 3) is irradiated with radiofrequency pulses, leading to an increase in the amplitude of the resonance of the first nucleus, as shown in Fig. 4. This effect indicates a dipolar interaction between the two spins (1 and 3). nOe is distance-dependent (1/rAX^6) and is crucial for structural biology applications.

nOe is essential for determining structural information in proteins and other biomolecules. It can be used to deduce distances between nuclei in the molecule, which helps in protein structure determination. It is particularly valuable for investigating the nature and register of secondary structures like ß-sheets.

The dependence of nOe on dynamics in a molecule is illustrated in Fig. 5. In proteins with varying flexibility, nOe can be qualitatively converted into distance information, but challenges like spin diffusion, where magnetization spreads along a chain of nuclei, need to be considered.

While nOe data can provide distance information, obtaining accurate distances can be challenging due to issues like spin diffusion. As a result, strategies have been established to prefer a larger set of qualitative nOes over a small set of precise distances, especially when interpreting NMR structures.

### Residual Dipolar Coupling (RDC)

In solution-state NMR, the dipolar coupling is typically averaged to zero due to the isotropic tumbling of molecules, which can make it challenging to detect. However, the presence of specific alignment media can induce a degree of molecular alignment, which, in turn, generates weak residual dipolar couplings (RDC) in the NMR spectrum

RDC data yield orientational information within a global frame of reference. This information can be used to determine the relative orientations of different parts of the molecule. While J-coupling provides angular information in terms of dihedral angles, RDC offers a complementary set of information that pertains to the orientation of molecular bonds in space.

Measuring and interpreting RDC data can be challenging. Finding suitable alignment media for a specific protein may require experimentation. The accurate measurement of numerous RDC values is time-consuming and requires careful data interpretation.

### NMR resonance assignment

The process of NMR resonance assignment is a critical step in nuclear magnetic resonance (NMR) spectroscopy when aiming to derive detailed structural information at the atomic level.

NMR resonance assignment can help detect protein-protein interactions. When the spectrum of an unlabeled protein B is titrated with a 15N-labeled protein A, changes in the NMR spectrum of B can indicate the interaction. However, other biophysical methods may provide similar information more cost-effectively.

Predicting NMR chemical shifts for a protein's resonances is extremely challenging, even if the protein's 3D structure is known. Therefore, experimental assignment is necessary.
Resonance assignment is like solving a jigsaw puzzle, where each resonance corresponds to a puzzle piece, and the goal is to connect them correctly.

### Two Resonance Assignment Strategies:

There are two main strategies for NMR resonance assignment:

***The 1H-1H approach*** : Utilizes two-dimensional NMR experiments, such as COSY, TOCSY, and NOESY, to establish connections between resonances within and between adjacent residues. This approach is suitable for small, well-behaved proteins.

***The 3D triple resonance (1H-15N-13C) approach*** : Introduced in the 1990s, it relies solely on J-coupling and is used for larger and more complex proteins.
it is well-suited for proteins of larger molecular weight because it copes with broader resonance lines and minimizes accidental resonance overlaps.


### NMR Tools

***Chemical Shift Perturbation (CSP)*** : CSP is a widely used NMR method that involves monitoring changes in chemical shifts in the NMR spectra of one molecule (usually a protein) as the concentration of the binding partner increases. This approach can provide information about the interaction interface between the two molecules.

***Paramagnetic Relaxation Enhancement (PRE)*** : Paramagnetic probes, such as nitroxide radicals, Mn2+ chelates, or lanthanides, can be introduced site-specifically in one of the interacting partners. These paramagnetic species induce paramagnetic relaxation enhancements (PRE) and pseudo-contact shifts (PCS) in the NMR spectra of the other partner. PRE effects primarily cause line broadening, while PCS leads to changes in resonance frequencies. PRE and PCS can provide long-range structural restraints and information about the relative orientation of interacting partners.

***Intermolecular Nuclear Overhauser Effect (nOe)*** : Intermolecular NOEs are NMR interactions that occur between nuclei of different molecules when they are in close proximity. These interactions can provide insights into the structure and dynamics of the complex formed by two interacting molecules.

***H/D Exchange Rates*** : Monitoring the exchange rates of hydrogen and deuterium (H/D) in amide protons upon complex formation can provide information about the interaction interface and the dynamics of the complex.

***Residual Dipolar Couplings (RDC)*** : RDC measurements can provide information about the orientation of internuclear vectors between interacting partners, particularly when anisotropic media are used to partially align the complex.


### NMR-Derived Model vs. X-ray Crystallography

In X-ray crystallography, a crystal of the protein is exposed to X-rays, and the resulting diffraction pattern is used to create an electron-density map, which is used to build an atomic model of the protein. In contrast, NMR-derived 3D structures are not images of the real structure but rather models of the structure that are consistent with experimental data.



### Usage

Biological NMR has evolved to encompass various applications beyond protein structure determination. While the number of NMR-solved structures has not grown as rapidly as those solved by X-ray crystallography, NMR offers unique advantages. It can investigate the dynamics of molecules over a wide range of time scales, making it useful for studying slow and fast molecular motions. In drug discovery, NMR is valuable for mapping chemical shifts to understand how proteins interact with ligands, making it a powerful tool for screening and optimizing potential drug candidates. 



1. NMR Parameters and Nuclei: NMR spectroscopy is used to study nuclei with a net spin, and hydrogen (1H) is particularly well-suited for NMR because it has a spin of 1/2. In contrast, carbon (12C), nitrogen (14N), and oxygen (16O) are not easily visible using NMR due to their isotopes. However, it's possible to enhance NMR sensitivity for these nuclei by "isotope-labeling" with isotopes like 13C and 15N, which has become more cost-effective over time.

2. Magnetic Field and Anisotropic Interactions: NMR experiments are conducted in a strong static magnetic field (B0), typically several Tesla, aligned along the +z axis. This magnetic field makes the environment non-isotropic, and all interactions experienced by the nuclear spins depend on the orientation of the molecule with respect to this magnetic field. Mathematically, these anisotropic NMR interactions are described by second-rank tensors or 3x3 matrices.

3. Liquid-State NMR: In liquid-state NMR, the molecule being studied is freely rotating with a correlation time (τc) typically ranging from 1 to 50 nanoseconds. This time is much smaller than the acquisition time for NMR data. When the rotation is isotropic, all interactions average out, and only the isotropic component of the interactions is observed. This property is why solution NMR spectra often display sharp resonances compared to solid-state NMR spectra, where the molecules are less mobile and the interactions are less averaged out.



### Conclusion


1. **Established Tool in Structural Biology:** NMR has become an established and powerful tool in the field of structural biology, offering atomic-level insights into the structure and dynamics of biological molecules, particularly proteins.

2. **Recent Technological Advances:** Recent advancements in isotopic labeling, magnet technology, electronics, and spectroscopy have pushed the boundaries of biomolecular NMR, allowing for more precise and comprehensive structural studies.

3. **Spectral Parameters for Conformational Information:** NMR spectral parameters, including J-coupling, nuclear Overhauser effects (nOe), residual dipolar couplings (RDC), and chemical shifts, provide detailed conformational information at the atomic level, enabling the determination of protein structures.

4. **Molecular Interaction Studies:** NMR spectroscopy is also a valuable tool for mapping molecular interactions, whether with small cofactors, RNA fragments, or other proteins. Spectral parameters are used to gain insights into these interactions.

5. **Intrinsically Disordered Proteins (IDPs):** NMR has played a crucial role in studying intrinsically disordered proteins (IDPs), which are unlikely to crystallize. It provides insights into their structural deviations, local structural propensity, and transient long-range contacts.

6. **Protein Dynamics:** NMR spectroscopy is used to investigate the dynamics of proteins over a wide range of time scales. This includes real-time NMR, exchange spectroscopy, and relaxation measurements to understand how proteins change over time.

7. **Limitations:** One of the limitations of NMR is the size of the proteins it can effectively study, although advances have allowed for the study of larger proteins, with some reaching up to 1 MDa.

8. **Emerging Areas:** Solid-state NMR and in-cell NMR spectroscopy are emerging fields in NMR research. Solid-state NMR is becoming a complementary tool for studying large soluble multimeric proteins and membrane-bound macromolecules. In-cell NMR is used to investigate the influence of complex and crowded cellular environments on protein structure and function.

9. **Complementary Perspective:** While NMR may not compete in terms of the number of proteins studied compared to high-throughput methods in proteomics, it offers an alternate and valuable perspective on many systems of biological interest.







