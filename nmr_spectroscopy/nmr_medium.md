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

![fig_2](https://github.com/LoqmanSamani/protein_structure_analysis/blob/systembiology/images/nmr_medium2.png)


To illustrate this concept, consider an Ethane molecule (refer to Figure 3). The carbon atoms in Ethane are equivalent in their arrangement, as are the hydrogen atoms. However, only the protons (hydrogen atoms) can contribute to an NMR signal. This is attributed to the fact that carbon-12, present in Ethane, has a Spin Quantum Number (I = 0), rendering it unable to produce an NMR signal. In contrast, hydrogen atoms possess a Spin Quantum Number (I) of -1/2 or +1/2, and only nuclei with I = 1/2 (either negative or positive) can generate an NMR signal. Consequently, as all protons in Ethane are equivalent, they collectively yield a single peak in NMR spectroscopy when subjected to an external magnetic field.

In the absence of an external magnetic field, nuclei have no discernible spin (spin is zero). However, when a magnetic field is applied, nuclei can adopt two spin states. The lower energy state is considered the ground state, while the higher energy state is designated the excited state. When electromagnetic radiation is supplied as energy, it promotes nuclei from the ground state to the excited state, resulting in observable signals in NMR spectroscopy.

![fig_3](https://github.com/LoqmanSamani/protein_structure_analysis/blob/systembiology/images/nmr_medium3.png)

#### Chemical Shift & Shielding Effect

Revisiting our illustrative molecule, Ethane (refer to Figure 3), let's scrutinize the behavior of the hydrogen atom at position 6 under the influence of an external magnetic field. This hydrogen nucleus, with a spin quantum number (I) equal to half, is capable of generating an NMR signal. However, its proximity to surrounding electrons introduces a phenomenon known as the ***shielding effect***.

Imagine applying a magnetic field (ß0) in a specific direction. The surrounding electrons, being negatively charged, respond by generating a magnetic field in the opposite direction, denoted as ***sigma * ẞ0 (δß0)***. The combined effect of these opposing magnetic fields results in an effective magnetic field ***ß(eff) = ß0 (1 - δ)***. In essence, the external magnetic field is not entirely transferred to the protons due to the shielding effect exerted by the surrounding electrons.

In NMR spectroscopy, the focus is not on measuring the shielding effect directly. Instead, it involves comparing this shielding effect with a standard, often ***TMS (tetramethylsilane)*** as depicted in Figure 4. The disparity between the standard (TMS) and the observed peak in the sample constitutes the ***chemical shift***. In NMR spectra, the shielding effect demonstrates an inverse relationship with the chemical shift – a higher shielding effect correlates with a lower chemical shift (figure 4).


![fig_4](https://github.com/LoqmanSamani/protein_structure_analysis/blob/systembiology/images/nmr_medium4.png)


#### Magnetic Anisotropy

In Figure 5, we examine a propane molecule, revealing distinct NMR signals for its three different protons (I, II, and III). Notably, protons I and III exhibit identical signals at 0.9 ppm, while proton II manifests a different signal at 1.3 ppm. Transitioning to another molecule, ethyl bromide (Figure 5), we again observe two distinctive signals for its proton types (I and II) in a proton NMR spectrum. However, the locations of these signals differ, with proton I at 1.7 ppm and proton II at 3.5 ppm.

An intriguing observation arises when comparing proton II in propane to proton II in ethyl bromide, where despite having an equivalent number of protons, their NMR signals diverge (1.3 ppm in propane and 3.5 ppm in ethyl bromide). This disparity can be attributed to the variance in electronegativity. Specifically, in ethyl bromide, carbon II is connected to a highly electronegative bromide atom, whereas in propane, carbon III is bonded to a hydrogen atom.

The high electronegativity of bromide induces a higher chemical shift for the second proton in ethyl bromide compared to its counterpart in propane. In essence, the magnetic anisotropy, influenced by electronegativity differences, manifests as varying chemical shifts in NMR spectra.

To summarize, the divergent electronegativities in molecular environments result in distinct chemical shifts, offering a nuanced understanding of magnetic anisotropy and its impact on proton NMR signals in molecules like propane and ethyl bromide.

![fig_5](https://github.com/LoqmanSamani/protein_structure_analysis/blob/systembiology/images/nmr_medium5.png)


#### Sample Preparation in NMR Spectroscopy

Sample preparation stands as a pivotal step in NMR spectroscopy, determining the success and accuracy of the experiment. It involves dissolving the test substance, the target for detection, in a suitable solvent to create a solution that can be analyzed within the spectrometer to obtain component spectra. The selection of an appropriate solvent is a critical aspect of this process.

Several key properties define an ideal solvent for NMR experiments:

- ***Purity***: The solvent should be available in a pure form to ensure that it does not introduce any impurities that could impact the spectra.

- ***Availability***: It should be readily accessible to facilitate ease of experimentation.

- ***Stability***: The chosen solvent should be stable under the experimental conditions, ensuring consistent and reliable results.

- ***Solvent Capacity***: A good solvent should effectively dissolve the sample substance to create a homogenous solution for analysis.

- ***Non-Interference***: The solvent should not interfere with the study. Given that NMR spectroscopy primarily targets nuclei with a spin quantum number (I) of 1/2, the solvent should ideally lack nuclei with I = 1/2 to prevent interference.

One suitable choice for non-polar samples is ***carbon tetrachloride***, which lacks nuclei with I = 1/2. Additionally, modified solvents like ***deuterated water (D₂O)*** and ***deuterated chloroform*** are commonly used. These solvents have modified hydrogen atoms, featuring one proton and one neutron in their nuclei, ensuring they do not interfere with the study.

Moreover, certain solvents contain nuclei with I = 1/2 but have NMR signals in a different range than the sample, preventing signal overlap and interference.

In conclusion, the careful consideration of solvent properties is crucial for successful NMR spectroscopy, guaranteeing accurate and unobstructed analysis of the target substance.


#### Fourier Transform in NMR Spectroscopy

The ***![Fourier Transform (FT)](https://github.com/LoqmanSamani/protein_structure_analysis/blob/systembiology/nmr_spectroscopy/fourier_transform.md)*** step in NMR (Nuclear Magnetic Resonance) is a pivotal mathematical process integral to the extraction of valuable information from raw NMR signals. In an NMR experiment, a sample is subjected to a strong magnetic field and radiofrequency pulses, prompting nuclear spins to absorb energy. As these spins relax back to their equilibrium states, they emit signals known as ***Free Induction Decay (FID)***, initially collected in the time domain.

The FID, a complex exponential signal with both amplitude and phase information, is represented as a series of data points sampled over time. To unveil the underlying frequencies of nuclear spins within the sample, the Fourier Transform is applied to convert this ***time-domain data*** into the ***frequency domain***(figure 6). This mathematical operation decomposes the complex FID into a spectrum, where each peak corresponds to a specific frequency of nuclear precession within the sample.

The resulting frequency-domain spectrum provides crucial insights. Peaks in the spectrum are associated with different atomic nuclei in the sample, and their positions reveal information about chemical shifts, indicating the chemical environment of the nuclei. Additionally, the splitting patterns in the spectrum convey details about nuclear couplings and the presence of neighboring nuclei.

![fig](https://github.com/LoqmanSamani/protein_structure_analysis/blob/systembiology/images/ft4.png)



#### Precession, Relaxation, and Resonance

Delving into the realm of magnetic resonance involves a profound exploration of atomic nuclei and their intriguing behavior within a magnetic field. To comprehend this, we must grasp two fundamental principles governing atomic nuclei.

Firstly, akin to tiny compass needles in a magnetic field, atomic nuclei possess inherent magnetism, aligning themselves with the field. This alignment forms the basis of the first principle.

The second principle delves into nuclear magnetism, emphasizing that nuclei possess angular momentum, resulting in a phenomenon known as "precession." Regardless of their orientation, nuclei exhibit a relatively consistent precession frequency. However, detecting this precession is most challenging when nuclei are either in a natural low-energy position or in an excited state. The precession becomes highly visible when nuclei exist in any intermediate orientation between the excited (pointing up) and natural (pointing down) states.

To initiate the appearance of a nucleus's precession frequency, a magnetic force must be applied at the same frequency as the precession. This synchronized movement of torque, precisely matching the nucleus's precession frequency, is termed "resonance." Resonance becomes a central concept in nuclear magnetic resonance, as it reorients atomic nuclei from their natural equilibrium states within a magnetic field.

In essence, nuclear magnetic resonance involves the application of a torque orthogonal to the magnetic field's torque on the spins. This torque must vary with time or oscillate at an exact frequency matching the natural precession frequency of nuclear spins. This synchronized application of torque, facilitated by a magnetic field, constitutes the phenomenon of nuclear magnetic resonance.

Regarding the element of time in this process, the gradual reorientation of nuclei occurs over time, impacting the visibility of precession. As time elapses, the precession's prominence diminishes, eventually leading to a return to some form of equilibrium.















### first order NMR spectra

it is the spectra where the chemical shift(Δv) is too greater than coupling constant(j) (Δv >> j), this type of spectra will be observed under the low resolution spectra condition(figure 5, a), and if we increase the resolution (high resolution spectra), the peak in the figure 5 a will be devided into two peak(it can be observed into two peak) (figure 5 b), the distance between this two spreaded peak is known as coupling constant(j) 

### chemical and magnetic equivalence



### rules of nmr
1) all magnetically equivalent nuclei will give one nmr signal
2) 










