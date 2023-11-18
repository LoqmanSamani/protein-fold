## Introduction to nuclear magnetic resonance (NMR)






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




#### Precession, Relaxation, and Resonance

Delving into the realm of magnetic resonance involves a profound exploration of atomic nuclei and their intriguing behavior within a magnetic field. To comprehend this, we must grasp two fundamental principles governing atomic nuclei.

Firstly, akin to tiny compass needles in a magnetic field, atomic nuclei possess inherent magnetism, aligning themselves with the field. This alignment forms the basis of the first principle.

The second principle delves into nuclear magnetism, emphasizing that nuclei possess angular momentum, resulting in a phenomenon known as ***precession***. Regardless of their orientation, nuclei exhibit a relatively consistent precession frequency. However, detecting this precession is most challenging when nuclei are either in a natural low-energy position or in an excited state. The precession becomes highly visible when nuclei exist in any intermediate orientation between the excited  and natural states.

To initiate the appearance of a nucleus's precession frequency, a magnetic force must be applied at the same frequency as the precession. This synchronized movement of torque, precisely matching the nucleus's precession frequency, is termed ***resonance***. Resonance becomes a central concept in nuclear magnetic resonance, as it reorients atomic nuclei from their natural equilibrium states within a magnetic field.

In essence, nuclear magnetic resonance involves the application of a torque orthogonal to the magnetic field's torque on the spins. This torque must vary with time or oscillate at an exact frequency matching the natural precession frequency of nuclear spins. This synchronized application of torque, facilitated by a magnetic field, constitutes the phenomenon of nuclear magnetic resonance.

Regarding the element of time in this process, the gradual reorientation of nuclei occurs over time, impacting the visibility of precession. As time elapses, the precession's prominence diminishes, eventually leading to a return to some form of equilibrium.



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

In the presence of an external magnetic field, nuclei exhibit distinct orientations. When subjected to a magnetic field along the Z-axis (Figure 2), nuclei can either rotate parallel to the applied magnetic field (lower energy and more stable, I = +1/2) or in the opposite direction (higher energy and less stable, I = -1/2). This phenomenon forms the basis of NMR spectroscopy, where the impact of the external magnetic field on these nuclear orientations is thoroughly examined.
In the realm of NMR, nuclei in the presence of a magnetic field exist in two states: a ground state and an excited state. When sufficient energy is supplied, nuclei in the ground state can transition to the excited state, resulting in absorption observed in the spectroscopy.


![fig_2](https://github.com/LoqmanSamani/protein_structure_analysis/blob/systembiology/images/nmr_medium_2.png)

<sub>**Figure 2:** Nuclei, when subjected to an external magnetic field (ß0), display distinct orientations. Nuclei with a spin quantum number (I) of +1/2 align parallel to the applied magnetic field, denoted as the aligned orientation, while those with I = -1/2 align in the opposite direction, forming the opposite orientation.</sub>


To illustrate this concept, consider an Ethane molecule (Figure 3). The carbon atoms in Ethane are equivalent in their arrangement, as are the hydrogen atoms. However, only the protons (hydrogen atoms) can contribute to an NMR signal. This is attributed to the fact that carbon-12, present in Ethane, has a Spin Quantum Number (I = 0), rendering it unable to produce an NMR signal. In contrast, hydrogen atoms possess a Spin Quantum Number (I) of -1/2 or +1/2, and only nuclei with I = 1/2 (either negative or positive) can generate an NMR signal. Consequently, as all protons in Ethane are equivalent, they collectively yield a single peak in NMR spectroscopy when subjected to an external magnetic field.

In the absence of an external magnetic field, nuclei have no discernible spin (spin is zero). However, when a magnetic field is applied, nuclei can adopt two spin states. The lower energy state is considered the ground state, while the higher energy state is designated the excited state. When electromagnetic radiation is supplied as energy, it promotes nuclei from the ground state to the excited state, resulting in observable signals in NMR spectroscopy.

![fig_3](https://github.com/LoqmanSamani/protein_structure_analysis/blob/systembiology/images/nmr_medium3.png)

<sub>**Figure 3:** The structural formula of Ethane. </sub> 


#### Chemical Shift & Shielding Effect

Revisiting our illustrative molecule, Ethane (Figure 3), let's scrutinize the behavior of the hydrogen atom at position 6 under the influence of an external magnetic field. This hydrogen nucleus, with a spin quantum number (I) equal to half, is capable of generating an NMR signal. However, its proximity to surrounding electrons introduces a phenomenon known as the ***shielding effect***.

Imagine applying a magnetic field (ß0) in a specific direction. The surrounding electrons, being negatively charged, respond by generating a magnetic field in the opposite direction, denoted as ***sigma * ẞ0 (δß0)***. The combined effect of these opposing magnetic fields results in an effective magnetic field ***ß(eff) = ß0 (1 - δ)***. In essence, the external magnetic field is not entirely transferred to the protons due to the shielding effect exerted by the surrounding electrons.

In NMR spectroscopy, the focus is not on measuring the shielding effect directly. Instead, it involves comparing this shielding effect with a standard, often ***TMS (tetramethylsilane)*** as depicted in Figure 4. The disparity between the standard (TMS) and the observed peak in the sample constitutes the ***chemical shift***. In NMR spectra, the shielding effect demonstrates an inverse relationship with the chemical shift – a higher shielding effect correlates with a lower chemical shift (figure 4).


![fig_4](https://github.com/LoqmanSamani/protein_structure_analysis/blob/systembiology/images/nmr_medium4.png)

<sub>**Figure 4:** Chemical shift of a sample relative to tetramethylsilane (TMS) as the standard, expressed in parts per million (ppm). The chemical shift values provide insight into specific regions in the spectrum associated with different functional groups or types of chemical environments. Peaks to the left (up-field) indicate lower chemical shift values, while peaks to the right (down-field) indicate higher values.</sub> 


#### J-coupling

***J-coupling or scalar coupling*** in NMR spectroscopy provides crucial structural information in one-dimensional spectra. This phenomenon results from the interaction of spin states through chemical bonds, causing the splitting of NMR signals. For protons, adjacent nuclei pointing towards or against the magnetic field lead to two signals per proton, revealing connectivity patterns.
Coupling to n equivalent (spin ½) nuclei generates an n+1 multiplet with intensity ratios following ***Pascal's triangle***. Additional spins cause further splittings, such as a doublet of doublets (dd) from coupling to two different spin ½ nuclei. Coupling between chemically equivalent or distant nuclei typically has no observable effect. Long-range couplings in cyclic and aromatic compounds produce complex splitting patterns.

In the proton spectrum of ethanol (figure 5), the CH3 group forms a triplet split by neighboring CH2 protons, while the CH2 is a quartet split by neighboring CH3 protons. Coupling with spin-1/2 nuclei like phosphorus-31 or fluorine-19 follows similar principles, providing insight into chemical environments and the number of neighboring NMR-active nuclei.
Coupling, combined with chemical shift and integration, reveals details about both chemical environment and neighboring nuclei. In complex spectra or non-hydrogen nuclei, coupling becomes essential for distinguishing between different nuclei.

![fig_5](https://github.com/LoqmanSamani/protein_structure_analysis/blob/systembiology/images/nmr_medium_5.png)

<sub>**Figure 5:** Variation in chemical shift for hydrogen atoms in ethanol attributed to neighboring chemical groups.</sub>

#### Magnetic Anisotropy

In Figure 6, we examine a propane molecule, revealing distinct NMR signals for its three different protons (I, II, and III). Notably, protons I and III exhibit identical signals at 0.9 ppm, while proton II manifests a different signal at 1.3 ppm. Transitioning to another molecule, ethyl bromide (Figure 6), we again observe two distinctive signals for its proton types (I and II) in a proton NMR spectrum. However, the locations of these signals differ, with proton I at 1.7 ppm and proton II at 3.5 ppm.

An intriguing observation arises when comparing proton II in propane to proton II in ethyl bromide, where despite having an equivalent number of protons, their NMR signals diverge (1.3 ppm in propane and 3.5 ppm in ethyl bromide). This disparity can be attributed to the variance in electronegativity. Specifically, in ethyl bromide, carbon II is connected to a highly electronegative bromide atom, whereas in propane, carbon III is bonded to a hydrogen atom.

The high electronegativity of bromide induces a higher chemical shift for the second proton in ethyl bromide compared to its counterpart in propane. In essence, the magnetic anisotropy, influenced by electronegativity differences, manifests as varying chemical shifts in NMR spectra.

To summarize, the divergent electronegativities in molecular environments result in distinct chemical shifts, offering a nuanced understanding of magnetic anisotropy and its impact on proton NMR signals in molecules like propane and ethyl bromide.

![fig_6](https://github.com/LoqmanSamani/protein_structure_analysis/blob/systembiology/images/nmr_medium_6.png)

<sub>**Figure 6** The observed differences in chemical shift values between protons in propane and ethyl bromide highlight the impact of magnetic anisotropy. </sub>


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

The FID, a complex exponential signal with both amplitude and phase information, is represented as a series of data points sampled over time. To unveil the underlying frequencies of nuclear spins within the sample, the Fourier Transform is applied to convert this ***time-domain data*** into the ***frequency domain***(figure 7). This mathematical operation decomposes the complex FID into a spectrum, where each peak corresponds to a specific frequency of nuclear precession within the sample.

The resulting frequency-domain spectrum provides crucial insights. Peaks in the spectrum are associated with different atomic nuclei in the sample, and their positions reveal information about chemical shifts, indicating the chemical environment of the nuclei. Additionally, the splitting patterns in the spectrum convey details about nuclear couplings and the presence of neighboring nuclei.



#### Mathematics behind Fourier Transform


***Representation of Intensity in the Frequency Domain***

The intensity at a specific frequency (f) is represented as the area under the product of Free Induction Decay (FID) and a cosine function:

    S_spectrum​(f) = ∫ S_FID​(t) ⋅ cos(2πft) dt

Here, S_FID(t) is the FID, and the integral is computed from 0 to infinity. However, in practice, the integration is performed up to a point where the FID becomes negligible.

***Discrete Fourier Transform in the Spectrometer***

In the spectrometer, the FID is represented by a set of sampled points. The Fourier Transform is applied by multiplying each FID point by the value of a trial cosine wave corresponding to that point in time. The resulting product functions are then summed to generate the intensity in the spectrum at the frequency of the trial cosine wave:

    S_spectrum​(f) = ∑ S_FID​(ti​) ⋅ cos (2πfti​)


Here, ti represents the time corresponding to the ith data point, S_FID(ti) is the value of the FID at that time, and NN is the total number of data points.

***Representation of the FID***

The x- and y-components of the magnetization generated by a 90° pulse about the x-axis can be represented as:

    My​ = −M0​ ⋅ cos(Ωt)

    Mx​ = M0​ ⋅ sin(Ωt)

If a 90° pulse about y were used, the components would switch:

    Mx​ = M0​ ⋅ cos(Ωt)

    My​ = M0​ ⋅ sin(Ωt)

The complex signal, S(t), representing both x- and y-components, is expressed as:

    S(t) = S0​ ⋅ exp(iΩt) ⋅ exp(−t/T2​)

This complex time-domain signal is a sum of terms if multiple resonances are present.

    S0 * exp(iΩt) * exp(-t/T2) + S0 * exp(iΩt) * exp(-t/T2) + ...

where each resonance i has its own intensity(S0), frequency(omega), and decay constant(T2)


***Fourier Transformation and Frequency-Domain Signal***

The Fourier transformation of the complex time-domain signal results in a complex frequency-domain spectrum:

    S(t) --(FT)--> ​S(ω) = S0​ ⋅ [A(ω)+i⋅D(ω)]

Here, A(ω) and D(ω) are Lorentzian absorption and dispersion mode lineshape functions, respectively.

***Arbitrary Phase in Time-Domain Signal***

The time-domain signal measured by the spectrometer has an arbitrary phase (Φ). The final representation of the frequency-domain signal includes this phase:

    S(t) = S0​ ⋅ exp(iΩt) ⋅ exp(−Rt) ⋅ exp(iΦ)

This results in both real and imaginary components in the frequency-domain spectrum.

#### Implementation of Fourier Transform using Python

```python

import numpy as np
import matplotlib.pyplot as plt



# Sample NMR time-domain data (FID - Free Induction Decay)

time = np.linspace(0, 1, 1000)  # Time values from 0 to 1 s
fid = np.sin(2 * np.pi * 10 * time) + 0.5 * np.sin(2 * np.pi * 20 * time)


# Perform Fourier Transform

spectrum = np.fft.fft(fid)
spectrum = np.fft.fftshift(spectrum)  # Shift the zero frequency component to the center



# Calculate the corresponding frequencies

sample_rate = 1.0 / (time[1] - time[0])
freq = np.fft.fftfreq(len(fid), 1 / sample_rate)
freq = np.fft.fftshift(freq)



# Plot the results

plt.figure(figsize=(10, 5))
plt.subplot(2, 1, 1)
plt.plot(time, fid)
plt.title("NMR Time-Domain Data (FID)")
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")

plt.subplot(2, 1, 2)
plt.plot(freq, np.abs(spectrum))
plt.title("NMR Frequency-Domain Spectrum")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Amplitude")


plt.tight_layout()
plt.show()

```


![fig](https://github.com/LoqmanSamani/protein_structure_analysis/blob/systembiology/images/ft4.png)

<sub>**Figure 7:** Fourier Transform Implementation in Python</sub>


### One-dimensional nuclear magnetic resonance spectroscopy (1D NMR)

In a one-dimensional NMR experiment, a sample is subjected to a strong external magnetic field(ß0), causing the nuclear spins of certain isotopes, commonly hydrogen (protons), to align either parallel or antiparallel to this field. ***Radiofrequency pulses*** are then applied, perturbing this alignment and prompting the spins to precess around the magnetic field. As these spins return to their equilibrium state, they emit signals, known as ***Free Induction Decay (FID)***, which are detected and translated into a spectrum.
The resulting one-dimensional NMR spectrum is a plot of signal intensity against the resonance frequency, known as the chemical shift. Each peak in the spectrum corresponds to a distinct set of protons in the molecule, providing information about their chemical environment. The chemical shift, measured in parts per million (ppm), reveals the magnetic environment of the protons and aids in identifying functional groups and molecular structures.
Additionally, the splitting patterns observed in the spectrum arise from spin-spin coupling, also known as J-coupling, elucidating the connectivity between neighboring protons. J-coupling patterns convey valuable information about the number of adjacent protons and their relative positions in the molecule.


### Two-dimensional nuclear magnetic resonance spectroscopy (2D NMR)

2D NMR constitutes a collection of techniques within nuclear magnetic resonance spectroscopy (NMR) that presents data plotted in a space defined by two frequency axes, as opposed to the single axis in one-dimensional NMR. Prominent methods within 2D NMR include ***correlation spectroscopy (COSY)***, ***J-spectroscopy***, ***exchange spectroscopy (EXSY)***, and ***nuclear Overhauser effect spectroscopy (NOESY)***. These multidimensional spectra offer enhanced insights into molecular structures, particularly for complex compounds challenging to analyze using one-dimensional NMR.

The standard framework for a 2D NMR experiment involves a sequence of radio frequency (RF) pulses interspersed with delay periods. These pulses' timing, frequencies, and intensities distinguish various NMR experiments. Typically, a 2D experiment encompasses four stages: the preparation period, creating a magnetization coherence through RF pulses; the evolution period, where nuclear spins precess freely without pulses; the mixing period, manipulating coherence into a state yielding an observable signal; and the detection period, where the free induction decay signal is observed, akin to one-dimensional FT-NMR.

The two dimensions in a 2D NMR experiment correspond to frequency axes representing a chemical shift. Each axis relates to one of the two time variables: the evolution period's duration, ***evolution time***, and the time elapsed during the detection period ,***detection time***. These time series are converted into frequency series through a ***two-dimensional Fourier transform***. A single 2D experiment is generated as a series of one-dimensional experiments, each with a distinct evolution time. The entire detection period is recorded in each experiment.

The outcome is a plot illustrating intensity values for each pair of frequency variables. Peak intensities are often represented using [***contour lines***](https://en.wikipedia.org/wiki/Contour_line) or distinct colors, providing a comprehensive representation of the molecular structure in 2D NMR spectra.



![fig_8](https://github.com/LoqmanSamani/protein_structure_analysis/blob/systembiology/images/nmr_medium_8.png)

<sub>**Figure 8:** COSY spectrum of Ethyl-2-butenoate</sub>





### NMR Spectroscopy in Protein Structure Determination

Nuclear Magnetic Resonance (NMR) spectroscopy has experienced significant advancements, particularly in the area of ***protein structure*** determination, marking an efficient technique in structural biology. The main goal of NMR studies of proteins is to obtain  ***high-resolution three-dimensional*** structures, similar to those  that can be achieved through X-ray crystallography. Despite having a limitation for proteins larger than 35 kDa, a boundary that has been surpassed in specific cases, NMR spectroscopy remains unparalleled for obtaining detailed structural insights into ***intrinsically unstructured proteins***, which are challenging for other methods.

In contrast to X-ray crystallography, NMR spectroscopy excels in capturing the dynamics and structural changes in proteins, making it a valuable tool for studying ***Conformation Activity Relationships (CAR)***. It plays an important role in comparing protein structures before and after interactions, such as with drug candidates, providing important information for ***drug development***.

The greater complexity resulting from the larger number of atoms in protein molecules compared to small organic compounds poses a challenge for spectral analysis. Basic ***1D spectra*** become convoluted due to ***overlapping signals***, making direct analysis impractical. To address this issue, ***multidimensional (2D, 3D, or 4D)*** experiments have been developed. ***Isotopic labeling***, particularly with ***13C and 15N***, becomes essential as it enhances the spectral resolution. The naturally occurring 12C isotope is not NMR-active, and the nuclear quadrupole moment of 14N impedes high-resolution information retrieval.

A fundamental approach in protein structure determination involves ***Nuclear Overhauser Effect (NOE)*** experiments. These experiments measure distances between atoms within a molecule, which provide important spatial information. Subsequently, the distances obtained are employed to solve a ***distance geometry problem***, yielding a high-fidelity 3D structure of the protein. The application of isotopic labeling significantly enhances the effectiveness of NOE experiments.

Beyond structural determination, NMR spectroscopy is a powerful tool for investigating the dynamics and conformational flexibility of various protein regions. This capability adds another layer of understanding to the functional aspects of proteins, offering insights into their behavior under different conditions.


  

### Applications of NMR in biology

Nuclear magnetic resonance (NMR) spectroscopy has proven to be a revolutionary tool in the field of biology, providing insights into the intricate world of protein structures. Its applications go far beyond the limits of traditional structural biology methods, providing a dynamic and detailed perspective on the behavior and conformation of biomolecules.

***Investigation of protein structures***

One of the most important applications of NMR in biology is to decipher the three-dimensional structures of proteins. In contrast to other techniques such as X-ray crystallography, NMR spectroscopy enables the investigation of proteins in their native environment in the solution phase. This is particularly important when studying proteins that are partially or completely unstructured and provides a comprehensive insight into their dynamic nature.

***Multidimensional NMR***

Innovation in protein NMR spectroscopy is highlighted by the development of multidimensional experiments. These include 2D, 3D and even 4D experiments designed to overcome the challenges posed by the complexity of protein structures. By isotopically labeling proteins with 13C and 15N, multidimensional NMR experiments enable the measurement of distances between atoms within a molecule and provide high-resolution structural information.

***Conformation-activity relationships (CAR)***

NMR spectroscopy plays a central role in the elucidation of conformation-activity relationships (CAR), particularly in drug discovery. It serves as a key technique for comparing the structures of proteins before and after interaction with drug candidates and provides insight into the structural changes underlying biochemical activities.

***Dynamics and flexibility***

Beyond static structures, NMR captures the dynamic nature of proteins. It provides information about the conformational flexibility and dynamics of different protein regions. This insight into the movements of proteins is crucial for understanding their biological functions, ligand binding and allosteric regulation.

***Groundbreaking studies***

Numerous  studies have highlighted the importance of NMR spectroscopy for biology and medicine. Examples include determining the structure of complex proteins that play a role in disease and revealing insights that serve as the basis for drug development. Notable contributions include the elucidation of the structure of membrane proteins, which pose a challenge to traditional methods of structural biology.

 
### Conclusion

In conclusion, Nuclear Magnetic Resonance (NMR) spectroscopy stands as a cornerstone in structural biology, providing unparalleled insights into the three-dimensional structures, dynamics, and interactions of biological macromolecules. Over its 80-year evolution, NMR has experienced transformative milestones, from the foundational principles of chemical shift and spin-spin interactions to the advent of multidimensional experiments and the elucidation of complex protein structures.

The key advantages of NMR spectroscopy lie in its ability to study proteins in solution, making it particularly valuable for intrinsically disordered proteins that defy crystallization. The technique's versatility spans from the elucidation of structural details to the exploration of molecular interactions, dynamics, and conformational changes crucial for understanding biological functions.

The journey of NMR spectroscopy has been marked by key milestones, including Nobel Prizes awarded for its foundational discoveries. The introduction of superconducting magnets, Fourier-transformed NMR, and two-dimensional techniques has significantly enhanced sensitivity and resolution, expanding the scope of NMR applications.

Basic principles, such as precession, relaxation, and resonance, underlie the functionality of NMR, while concepts like chemical shift, J-coupling, and magnetic anisotropy contribute to the richness of information obtained. Sample preparation, with careful consideration of solvent properties, and the crucial Fourier Transform step are integral aspects of successful NMR experiments.

In the realm of protein structure determination, NMR spectroscopy has emerged as a powerful tool, excelling in capturing dynamics and conformational changes. From high-resolution structures to the investigation of conformation-activity relationships, NMR plays a central role in drug development and functional studies.

In conclusion, the evolution of NMR spectroscopy has shaped it into an indispensable technique in the toolkit of structural biologists, providing a holistic approach to unraveling the mysteries of biological macromolecules.


----------------------------------------------------------------------

----------------------------------------------------------------------



### References

- Understanding NMR Spectroscopy by James Keeler.     
 

- Marion D. An introduction to biological NMR spectroscopy. Mol Cell Proteomics. 2013 Nov;12(11):3006-25. doi: 10.1074/mcp.O113.030239. Epub 2013 Jul 6. PMID: 23831612; PMCID: PMC3820920.
   

- Curry S. Structural Biology: A Century-long Journey into an Unseen World. Interdiscip Sci Rev. 2015 Jul 3;40(3):308-328. doi: 10.1179/0308018815Z.000000000120. Epub 2015 Dec 7. PMID: 26740732; PMCID: PMC4697198.
     

- Campbell, Iain D. ‘The Evolution of Protein NMR’. 1 Jan. 2013 : 245 – 264.
     

- Foundations of Structural Biology By Leonard J. Banaszak.


- Three-dimensional triple-resonance NMR spectroscopy of isotopically enriched proteins Kay, Lewis E. ; Ikura, Mitsuhiko ; Tschudin, Rolf ; Bax, Ad.
   

- Wüthrich K. The way to NMR structures of proteins. Nat Struct Biol. 2001 Nov;8(11):923-5. doi: 10.1038/nsb1101-923. PMID: 11685234.


- Fuloria NK, Fuloria S (2013) Structural Elucidation of Small Organic Molecules by 1D, 2D and Multi Dimensional-Solution NMR Spectroscopy. J Anal Bioanal Tech S11:001. doi: 10.4172/2155-9872.S11-001.



------------------------------

***Loghman Samani***  

***November 18, 2023***