# Nuclear Magnetic Resonance (NMR)

#### NMR_Spectroscopy is an analytic technique used to identify the molecular structure of a molecule. This is accomplished by analysing the magnetic fields around atoms within a molecule. 

## Terminology


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

- I (1/2 or -1/2);

Then:

- M = 2 * 1/2 + 1 = 2;

Then:

- S(+1/2) = ħ / (2 * π * 1/2) = ħ / (4 * π) and S(-1/2) = -ħ / (4 * π);

Then:

- μ(+1/2) = λ * S = +λ * ħ / (4 * π) and μ(-1/2) = λ * S = -λ * ħ / (4 * π);

Then:

- E(+1/2) = μ(+1/2) * ß⁰ = +λ * ħ * ß⁰ / (4 * π) and E(-1/2) = μ(-1/2) * ß⁰ = -λ * ħ * ß⁰ / (4 * π);

Then:

- ΔE = E(-1/2) - E(+1/2) = -λ * ħ * ß⁰ / (2 * π) (the minus sign here does not matter; only the magnitude is important);

Then:

- ΔE = ħν = λ * ħ * ß⁰ / (2 * π);

Then:

- ν = λ * ß⁰ / (2 * π).

Atoms with I equal -1/2 or I equal +1/2 can be excited, and so they can release an NMR signal. The question here is which atoms can be excited? The answer depends on the surrounding of the atom. For example, if we have a molecule H-C-Cl; if we study the nuclei structure of carbon (assuming it is carbon-12), it contains 6 protons and 6 neutrons, so the I should be 0, and it will not release any signal. However, due to the higher electronegativity of Cl, it attracts the shared electrons between it and carbon, so the carbon might be excited and release an NMR signal.

## Some Examples

### Molecule: Ethane



---------------------

        H   H
        |   |
    H - C - C - H
        |   |
        H   H

---------------------

- A molecule like ethane has two types of atoms (hydrogen and carbon).
- Carbon atoms cannot give signals.
- Hydrogen can give an NMR signal.


### Molecule: (a), (b), (c)

---------------------------------

        (a)      (b)     (c)
         H       H       H
         |       |       |
     H - C - H - C - H - C - Cl
         |       |       |
         H       H       H

---------------------------------

- Hydrogens (a), (b), and (c) are attached to different groups, so they all have different NMR signal strengths.

In NMR spectroscopy, the first step (after sample preparation) is to apply an external magnetic field (ß⁰). The applied magnetic field excites hydrogen atoms. The difference between two states (excited and unexcited state) depends on the surrounding electrons in an atom, which can show a "shielding effect." This effect is related to the strength of the released signal. It means that not all the energy from the external magnetic field (ß⁰) can be applied to the protons in the nuclei of the atom (in our case, carbons cannot be excited, but hydrogen can). Some protons align in the same rotation direction as ß⁰ (called aligned orientated Hs with I = +1/2), while others align in the opposite direction (opposite orientated Hs with I = -1/2). The next step is to supply energy in the form of EMR that resonates with the energy gap.

In an NMR process, the difference between the shielding effect (higher electronegativity, higher shielding effect, lower the chemical shift) and the standard (TMS), which has the highest shielding effect, will be measured.

-----------------------------------------------------------------------------------------------

                                    
                                       
       down-field                                                        up-field
                                                                   
        (carboxylic)   (aldehyde)   (aromatic)   (halides)   (alkanes)     (TMS)                                                
        ________________________________________________________________________
        12         10          8            6           4           2          0
                                      ppm
   





    down-field                                                               up-field
                                                       sample            TMS
                                                         |       d        |
        _________________________________________________|________________|
        12               ...        6           4            2          0
                                      ppm

    d = chemical shift: the difference between TMS and Sample  

-----------------------------------------------------------------------------------------------





Loghman Samani

04.November.2023


