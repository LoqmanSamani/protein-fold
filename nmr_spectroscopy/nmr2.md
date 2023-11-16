### Precession, Relaxation, and Resonance

In thinking about magnetic resonance we will be thinking about atomic nuclei and their behaviour in a magnetic field.
The first principle of understanding atomic nuclei in a magnetic field is to realize that those nuclei have magnetism and like a little compass-needles which is placed in a magnetic field they tend to line up with the field.
The second principle of understanding atomic nuclei (nuclear magnetism) is the idea that nuclei themselves have angular momentum and this angular momentum of nuclei makes an effect which is called "precession". actually it does not matter what orientation a nucleus has, the precession frequency is pretty much the same however if the nucleus is in the natural low energy position (not excited) or is excited the precession will be very hard to detect.
but anywhere in between pointing up(excited nucleus) and pointing down (natural state or not excited) the precession is very visible.
To make the precession frequency of a nucleus to appear a force (magnetic force) should be applied (at the same frequency as the precession frequency), so the movement of the torque following exactly the precession frequency of the nucleus is known as resonance and the frequency is exactly the same in both cases. The idea of resonance is really central to what we will be looking at and will be reoriented atomic nuclei from their natural equilibrium states in a magnetic field.

in summary when we perform nuclear magnetic resonance we have to apply a torque that is orthogonal to the torque of the magnetic field on the spins and we have to have that torque varying with time or oscillating an exact frequency match with natural precession frequency of the nuclear spinsand we apply this torque with a magnetic field and that is nuclear magnetic resonance.
we talk about time but how exactly time affects the resonance?  with time gradually nuclei will reorient and the visibility of the precession starts to die away and this returning some sort of equilibrium.


### Quantum Behaviour of Atomic Nuclei

If we imagine a nucleus (hydrogen nucleus known as a proton) as a ball and a vector inside of it which represent the direction of the angular momentum and the magnetism. In a magnetic field the lowest energy state is with magnetic dipole of the nucleus aligned with the magnetic field direction. Hydrogen nucleus (proton) has a two possible quantum states in the presence of magnetic field, spin up (low energy state) and spin down (high energy state). In general being quantum mechanics it is possible to have a coherent superposition of this states, in which one nucleus iss both spin up and spin down simultneously.
That is superposition that applies when we are observing protons precessing in a magnetic resonance experiment.

![nmr1](https://github.com/LoqmanSamani/protein_structure_analysis/blob/systembiology/images/nmr1.png)

When we have lots of nuclei in thermal eqilibrium with a magnetic field. At finite temperature all nuclei are not in the lower energy state a proportion of them are in the high energy state depend on the temperature but none of these two position are superposition state, the ratio of the population of each state is given by Boltzmann Factor (B = e ^(- yß0h/kBT)) when we decrease the temperature the ratio of lower enerpy spins will increase and in the zero (kelvin) al the nuclei will be in the lower energy state.(nmr2)
another way to incease the ratio of lower energy spins is to increase the applied magnetic field (ß0). in this case the thermal energy remains fixed.

![nmr2](https://github.com/LoqmanSamani/protein_structure_analysis/blob/systembiology/images/nmr2.png)

lets take our thermal magnetization equilibrium and see what we can do using the nuclear magnetic trick! like i said before the precession is invisible if the magnetization pointed along the vertical axis and we need to reorient the nuclei. once we do that precession will be visible we call this precession Larmor Frequency (ω = γß0) this equation is the most important equation in NMR. it says the precession frequency is proportional to the field at that the constant of proportionality gamma (gyromagnetic ratio or the magnetogyric ratio) depends on the nature of the nucleus. Hydrogen nuclei have the biggest gamma among the most sensitive nuclei to used for NMR. finally note the use of a coil(nmr3) pick up the induced voltage from the precessing spins, this constitude our free induction decay signal in NMR.

![nmr3](https://github.com/LoqmanSamani/protein_structure_analysis/blob/systembiology/images/nmr3.png)

summany: the strength of the magnetic field determines the amount of magnetization that we have when the atomic nuclei are in thermal equilibrium. The strength of the magnetic field also determines the precession frequency of the atomic nuclei, the higher the field the higher the precession frequency. wwe need to have a coil which used to produce an ascillating magnetic field in resonance with the nuclei that causes that talk which tips the mechanization out of the equilibrium into the transverse plane (the coil is used not only to transmit to the nuclei but also to receive the signal to pick up that electro-magnetic force that oscillating voltage as the nuclear magnetization processes around in the coil with signal gradually decaying away with time as the spins come back to the thermal equilibrium and we call that signal which decays with time "free induction decay").

![nmr4](https://github.com/LoqmanSamani/protein_structure_analysis/blob/systembiology/images/nmr4.png)



###  Acquiring a Free Induction Decay (FID)

Two principle things which we need for magnetic resonence: 1) ***A homogeneous static magnetic field***, it is provided by the earth magnetic field, which is about 60 micro Tesla and leads to allow more precission frequency for hydrogen nuclei just over of couple of kilohertz (it has to do with the sharpness of the resonance the degree of resolution that is possible with nuclear magnetic resonance) 2) ***An oscillating transverse magnetic field***, which is provided by a machine.
If we use water as a sample for an experiment (water has two hydrogen and one oxygen, the oxygen atom has actually no magnetizem and also no angular momentum and it means plays no part in the nuclear magnetic resonance phenomenon but the hydrogen(proton) around the oxygen atom have all we need) a sample of water has a gigantic numer of proton in it which are appropriate for the experiment.
there are actually lots of other molecules which can be use as a sample instead of water (molecules which contain hydrgon atoms in them such aliphatic chains which are in form of lipid molecules in our body. actually any molecules which contains hydrogen and can befound in liquid form is appropriate). first in this experiment we apply a prepolarizing pulse to produce some larger magnetization and then we have to disturb the spins from their thermal equilibrium state and for that we use the transverse oscillating magnetic field, we apply it for a short period of time as a pulse. what we detect is the excitation and the decay of that oscillating signal at the free induction decay over a period of time (in this case 2 seconds) (time domain curve) if we perform a fourier transform on that time domain data we get the spectrum and that spectrum lies in the frequency domain correcponds (in this case) to a definite peak.(nmr5)

![nmr5](https://github.com/LoqmanSamani/protein_structure_analysis/blob/systembiology/images/nmr5.png)


###  Field Homogeneity

***Why the homogeneity of the field is important?*** THe width of the peak in the spectrum related to the rate of decay of the signal in the frequency graph(FID), the faster the decay the broader the spectrum, so why the signal decays? initially the spins do have to go back to thermal equilibrium along the magnetic field so it has to be some die away of the signal due to this process, but in fact the rate of decay turns out to be a lot faster than that which would expect from return to thermal equilibrium. The reason of the fast decay is related to field homogeneity. Explanation: supposed we have a field that is perfectly homogeneous and the atomic nuclei are in different positions in that magnetic field and they will have exactly the same larmor precession frequency so the magnetization vectors associated with those atomic nuclei will precess at exactly the same rate, now imagine the magnetic field is inhomogeneous, that means it varies from place to place across the sample so atomic nuclei in different positions have different larmor precession frequencies, that means when we do the free induction decay experiment the magnetization vectors start together in phase but of the different precession frequencies the vectors gradually spread out so the total magnetization which is a sum of the different vectors is smaller as a result of them spreading out.


### Methods to prevent quick decay (improve homogeneity)

***Shimming*** : By putting a small amount of current in the magnetic coils (this method will be applied by the nmr-machine) we can correct for imperfections in the Earth's magnetic field.

***Spin Echo*** : oscillating transverse magnetic field has the effect of causing the nuclear magnetization to move from its equilibrium position pointing along the magnetic field to being precessing at some particular angle and the most favour angle maximum signal is that which have that magnetization lying in the transverse plane how do we get just the right amount of turning?  we have to apply the oscillating magnetic field for just the right amount of time, and that time produces that is known as a 90-degree pulse.(nmr6) for doing this we turn the magnetization vectors, which have different precession frequencies(some are faster and some slower), they will spread out with time but at some point we turn them over (180 degree), so they eventually come back into step with each other.

![nmr6](https://github.com/LoqmanSamani/protein_structure_analysis/blob/systembiology/images/nmr6.png)


### Magnetic Resonance Imaging (MRI)

For understanding the idea behind MRI and to be able to interpret the result of an MRI imagine we have a simple sample of water (this sample contains two tube of water which are placed parallel to each other). we said before that the precession frequency of nuclear spins depends on the strength of the magnetic field. imagine if we could make the magnetic field higher on the top tube and lower on the bottom tube, that would mean that the top tube would have a higher frequency and the bottom tube would have the lower frequency and if we look at the spectrum of that we should see two signals at two different frequencies corresponding to the two tubes of water, now imagine if we could make that magnetic field very linearly along the axis where the two tubes are separated, in that case we would have a very simple relationship between frequency and position in fact the frequency axis would also be a position axis, now the question is how could we do that? for doing it we use the ability of magnetic-coils inside the MRI machine, they could actually produce a variation of the magnetic field in all three orthogonal directions (along x, y and z axis), in this case we only need to have a variation of field along the z axis, so we are going to be able to separate those two tubes of water, so the idea is to use these gradient coils to produce the variation of magnetic field along the z axis in combination with a spin echo experiment(nmr7).

![nmr7](https://github.com/LoqmanSamani/protein_structure_analysis/blob/systembiology/images/nmr7.png)







