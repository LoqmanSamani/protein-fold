## Introduction into Nuclear Magnetic Resonance (NMR)  


<img src="https://github.com/LoqmanSamani/protein_structure_analysis/blob/systembiology/images/nmr-spectroscopy.jpg" width="1000">

*Figure 1: Picture of an NMR spectrometer*


### Introduction

In the 20th century, a variety of ***experimental techniques*** were used in an attempt to decipher the complicated three-dimensional structures of biological molecules. Among these, ***X-ray crystallography***, ***nuclear magnetic resonance (NMR)***, and ***electron microscopy*** stood out as pioneers in revolutionizing ***structural biology***. The discovery and application of ***X-rays*** to ***protein crystals*** marked a transformative era, enabling scientists to unveil the atomic details of biological molecules.

With the progress of these techniques, ***structural biology*** developed into an independent discipline that is intertwined with molecular biology, biochemistry and biophysics and is dedicated to deciphering the molecular architecture of biological macromolecules.

This article deals with the advantages of NMR spectroscopy and its central role in the elucidation of the ***three-dimensional structures of proteins***. NMR spectroscopy has proven to be a powerful tool that captivates biologists working in structural biology, providing unparalleled insights into the dynamics, interactions and structures of biological macromolecules at the atomic level. Beyond pure structural elucidation, NMR spectroscopy is proving invaluable in the study of molecular interactions by revealing crucial details about binding interfaces, kinetics and thermodynamics.

In particular, NMR plays a central role in the study of ***intrinsically disordered proteins (IDPs)***, a category that cannot be crystallized. The versatility of NMR spans various biological fields and demonstrates its capabilities in unraveling the mysteries of protein dynamics that are crucial for understanding their biological functions.


### Advances in NMR Spectroscopy

The inception of nuclear magnetic resonance (NMR) spectroscopy, approximately 80 years ago at the conclusion of World War II, marked a transformative moment in scientific exploration. Independently discovered by two distinguished groups led by Felix Bloch et al. and Edward Purcell et al., NMR's evolution has been punctuated by key milestones and innovations.

In the early years following its discovery, NMR found its way into the realm of chemistry when W.G. Proctor and F.C. Yu stumbled upon a serendipitous observation—distinct signals from the two nitrogens in NH4NO3. This marked the pioneering recognition of the ***chemical shift*** phenomenon, a cornerstone in NMR spectroscopy.

The narrative advanced with the detection of spin-spin interactions, known as the ***nuclear Overhauser effect***, by Overhauser in 1953. This breakthrough revelation showcased the communication between spins, observed during the saturation of electrons in metals.

Initially constrained by limited sensitivity, NMR faced a pivotal turning point with the advent of stronger electromagnetic fields. Superconducting magnets, introduced in the early 1960s, propelled NMR frequencies to 100 MHz for the 1H nucleus. Overcoming technical challenges related to ***magnet homogeneity*** and stability, the first 200 MHz spectrum of ethanol was triumphantly published in 1964.

A quantum leap in sensitivity followed with the introduction of ***Fourier-transformed (FT) NMR*** in 1966 by Ernst and Anderson. This methodological breakthrough allowed the simultaneous excitation and subsequent unravelling of all signals, paving the way for the development of myriad pulse sequences.

The journey into structural biology commenced in the 1960s with studies on small biological molecules like common amino acids. Synthetic and natural peptides, as well as paramagnetic proteins such as cytochrome c and myoglobin, dominated research endeavors during this period.

A pivotal moment arrived in 1976 with the introduction of ***two-dimensional (2D) NMR*** by R. Ernst et al., building upon the insightful idea of J. Jeener, a Belgian physicist.

The landscape shifted in 1985 with the publication of the first structure of a small ***globular protein***. Despite initial skepticism from the X-ray crystallography community, this marked the advent of three-dimensional (3D) NMR. This technique initially unfolded on unlabeled proteins in 1990, swiftly followed by a suite of triple resonance experiments using 15N and 13C labeled samples.

Today, NMR stands as a unique and versatile technique, offering a plethora of multidimensional experiments. 

#### Key Milestones in the Evolution of NMR Spectroscopy

- ***1952*** Nobel Prize in Physics awarded to Felix Bloch and Edward Purcell for their independent discovery of NMR.

- ***1950*** Proctor and Yu make the first observation of the ***chemical shift***.

- ***1952*** Production of the first commercial NMR spectrometer.

- ***1953*** Overhauser observes the ***nuclear Overhauser effect***.

- ***1960s*** Introduction of ***superconducting magnets***.

- ***1966*** Introduction of ***Fourier-transformed (FT) NMR***.

- ***1976*** Introduction of ***two-dimensional (2D) NMR*** by R. Ernst and others.

- ***1985*** Publication of the first structure of a small globular protein.

- ***1991*** Nobel Prize in Chemistry awarded to Richard Ernst for contributions to ***high-resolution NMR spectroscopy***.

- ***2002*** Nobel Prize in Chemistry awarded to Kurt Wüthrich for advancements in understanding the three-dimensional structures of biological macromolecules in solution.



Basic Principles of NMR Spectroscopy:

    Explain the fundamental principles of nuclear magnetic resonance.
    Discuss the interaction of magnetic nuclei with external magnetic fields.
    Introduce concepts like chemical shift, spin-spin coupling, and relaxation times.



### Basic Principles of NMR Spectroscopy


#### The Chemistry Behind Magnetic Resonance


- **Spin Quantum Number (I):** Possible values -1/2, 0, 1/2, 1 (atoms with I equal to 1/2 plus or minus can be excited).

- **External Magnetic Field Strength (ß⁰):**

- **EMR Energy:** In most cases, it uses radio-waves.

- **Magnetic Momentum Number (M):** The number of spin states of a nucleus in a magnetic field (M = 2I + 1).

- **Gyromagnetic Ratio (λ):**

- **Spin Angular Momentum (S):** S = ħ / (2 * π * I).

- **Magnetic Momentum (μ):** μ = λ * S.

- **Magnetic Energy (E):** E = μ * ß⁰.

- **Larmor Frequency or Larmor Precession (ν):** ν = λ * ß⁰ / (2 * π).


If:

    I = 1/2 

    or

    I = -1/2

Then:

    M = 2 * 1/2 + 1 = 2

Then:

    S(+1/2) = ħ / (2 * π * 1/2) = ħ / (4 * π)

    S(-1/2) = ħ / (2 * π * -1/2) = -ħ / (4 * π)

Then:

    μ(+1/2) = λ * S = +λ * ħ / (4 * π) 

    μ(-1/2) = λ * S = -λ * ħ / (4 * π)

Then:

    E(+1/2) = μ(+1/2) * ß⁰ = +λ * ħ * ß⁰ / (4 * π)  # The lower energy state (ground state)
   
    E(-1/2) = μ(-1/2) * ß⁰ = -λ * ħ * ß⁰ / (4 * π)  # The higher energy state (excited state)

Then:

    ΔE = E(-1/2) - E(+1/2) = -λ * ħ * ß⁰ / (2 * π) 
   
    The minus sign here does not matter, only the magnitude is important.

Then:

    ΔE = ħν = λ * ħ * ß⁰ / (2 * π)  # The energy gap between the two energy states

Then:

    ν = λ * ß⁰ / (2 * π)  # Larmor Frequency or Larmor Precession



#### How NMR Works: Understanding Spin States and Energy Transitions

In the presence of an external magnetic field, nuclei exhibit distinct orientations. When subjected to a magnetic field along the Z-axis (as illustrated in Figure 2), nuclei can either rotate parallel to the applied magnetic field (lower energy and more stable, I = +1/2) or in the opposite direction (higher energy and less stable, I = -1/2). This phenomenon forms the basis of NMR spectroscopy, where the impact of the external magnetic field on these nuclear orientations is thoroughly examined.
In the realm of NMR, nuclei in the presence of a magnetic field exist in two states: a ground state and an excited state. When sufficient energy is supplied, nuclei in the ground state can transition to the excited state, resulting in absorption observed in the spectroscopy.

![fig_2]()

To illustrate this concept, consider an Ethane molecule (refer to Figure 3). The carbon atoms in Ethane are equivalent in their arrangement, as are the hydrogen atoms. However, only the protons (hydrogen atoms) can contribute to an NMR signal. This is attributed to the fact that carbon-12, present in Ethane, has a Spin Quantum Number (I = 0), rendering it unable to produce an NMR signal. In contrast, hydrogen atoms possess a Spin Quantum Number (I) of -1/2 or +1/2, and only nuclei with I = 1/2 (either negative or positive) can generate an NMR signal. Consequently, as all protons in Ethane are equivalent, they collectively yield a single peak in NMR spectroscopy when subjected to an external magnetic field.

![fig_3]()

In the absence of an external magnetic field, nuclei have no discernible spin (spin is zero). However, when a magnetic field is applied, nuclei can adopt two spin states. The lower energy state is considered the ground state, while the higher energy state is designated the excited state. When electromagnetic radiation is supplied as energy, it promotes nuclei from the ground state to the excited state, resulting in observable signals in NMR spectroscopy.


#### Chemical Shift & Shielding Effect

Revisiting our illustrative molecule, Ethane (refer to Figure 3), let's scrutinize the behavior of the hydrogen atom at position 6 under the influence of an external magnetic field. This hydrogen nucleus, with a spin quantum number (I) equal to half, is capable of generating an NMR signal. However, its proximity to surrounding electrons introduces a phenomenon known as the ***shielding effect***.

Imagine applying a magnetic field (ß0) in a specific direction. The surrounding electrons, being negatively charged, respond by generating a magnetic field in the opposite direction, denoted as ***sigma * ẞ0 (δß0)***. The combined effect of these opposing magnetic fields results in an effective magnetic field ***ß(eff) = ß0 (1 - δ)***. In essence, the external magnetic field is not entirely transferred to the protons due to the shielding effect exerted by the surrounding electrons.

In NMR spectroscopy, the focus is not on measuring the shielding effect directly. Instead, it involves comparing this shielding effect with a standard, often ***TMS (tetramethylsilane)*** as depicted in Figure 4. The disparity between the standard (TMS) and the observed peak in the sample constitutes the ***chemical shift***. In NMR spectra, the shielding effect demonstrates an inverse relationship with the chemical shift – a higher shielding effect correlates with a lower chemical shift (figure 4).


![fig_4]()

















### first order NMR spectra

it is the spectra where the chemical shift(Δv) is too greater than coupling constant(j) (Δv >> j), this type of spectra will be observed under the low resolution spectra condition(figure 5, a), and if we increase the resolution (high resolution spectra), the peak in the figure 5 a will be devided into two peak(it can be observed into two peak) (figure 5 b), the distance between this two spreaded peak is known as coupling constant(j) 

### chemical and magnetic equivalence



### rules of nmr
1) all magnetically equivalent nuclei will give one nmr signal
2) 










