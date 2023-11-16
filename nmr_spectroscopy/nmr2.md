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


### Introduction to k-space

in this section we speak about ideas around spatial frequency which are important in understanding magnetic resonance imaging. we begin with the magnetic field ß0 (fig. nmr8, a) and the magnetization vector for our topic nuclei in thermal equilibrium. to start the magnetic resonance experiment is to tip that magnetization vector away for equilibrium which we do with the oscillating transverse magnetic field. as the spins precess around (fig. nmr8, b) they have different orientation in the transverse plane and in physics the way we describe that is the idea of phase, phase tell us the angle of the magnetization vector has in that precession in the transverse plane. the mathematical description of phase is something called phase vector (exp(iγß0t)), it is an exponential function which in the argument has the product of the larmor frequency γ ß0 and the time (t) and the product between frequency and time makes an angle in radians so the phase angle is what appears in the argument of the exponential function (exp(iγß0t)). 
in an MRI we have a magnetic field which varies with position, this system produce a variation of magnetic field along the vertical axis where there is a linear depedance on position. so we rewrite the larmor equation (ω = γß0) to describe the frequency of precession of the spins in this case (ω = γß + γG) the aditional term (γG) is because of the additional gradient coil field and this term adds frequency depends upon position along the z axis and it tells us in which the magnetic field varies along the z axis, so the result of this equation will be an angular frequency in radians per second. the atomic nuclei (the sample) are in different positions along the axis of the gradient coil, a set of magnetization vectors which all start together with the same phase but they distributed along that the axis at different positions (fig. mnr8, c), so what that will mean is that each of those magnetization vectors is going to experience an additional filed from the gradinet coil which will vary according to where they are, and so the subsequent precession will accure at different rates. the result will be a helix of phase in the spin system.

![nmr8](https://github.com/LoqmanSamani/protein_structure_analysis/blob/systembiology/images/nmr8.jpg)

In the picture (fig. nmr9) we can see that after the 90-degree pulse of the oscillating transverse magnetic field, the magnetization vectors of all atomic nuclei are tipped into the transverse plane, subsequently with time they will precess all of the same rimel frequency and we will have the free induction decay which will have a reducing episode with time because of relaxation. suppose we have a magnetic field gradient that we switch on at some period of time after the 90-degree pulse, because of this magnetic field the nuclear magnetization vectors all start off together in phase all in step with each other and because of the presence of the magnetic field gradient they start to acquire different phases depending on where they are and those different phases cause a helical variation and phase along the vertical axis because the magnetic field gradient is a field that is proportional to position and as time goes on the helix whines tighter and the pitch of the helix gets shorter. we can describe the change in wave lengths mathematically (K = γGt and  K = 2 pi / γ). If we obtain a signal at some perticular point in time, the signal would come about from the sum of the contributions of all the atomic nuclei, so the signal would have involved a summation which represented by this integral  ( S(k) = ∫ p(z) * exp(ikz) dz ), in this formula we summ up the contribution of each individual spin at different positions. they may not be equal number of atomic nuclei difference position and that would be the reason that we carry out in the imaging experiment to try to find out the distribution of an atomic nuclei in space.

![nmr9](https://github.com/LoqmanSamani/protein_structure_analysis/blob/systembiology/images/nmr9.png)

when we look at one dimensional imaging (fig.10) what we see is the signal (S(k)) we acquired in time also corresponds to a signal in k space because k is proportional to time and when we fourier transform ( S(k) = ∫ p(z) * exp(ikz) dz , here the transformation converts the acquired signal into density in real space) that to get the frequency domain we obtain the image.

![nmr10](https://github.com/LoqmanSamani/protein_structure_analysis/blob/systembiology/images/nmr10.png)


### k-space in multiple dimensions

Here we extend the idea of reciprocal space or k-space into two or more three dimensions, for doing that we need more than one gradient coil, which were represented in the one dimensional k-space, that means the gradient can be represented as a vector (G = ∇ß) depending on the direction of the gradient applied  by each individual coil and the Larmor precession equation in this case will be represented as (ω = γß + γG . r, here we have a dot product between gradient vector(γG ) and the position vector(r))(fig. nmr11).

![nmr11]()

In this case if we look at the pulse sequence, that we would use, the normal time course of events starting with a 90-degree pulse and seeing how thing evolve with time when we are using two- or three-dimensional gradient(fig. nmr12). Here we rewrite the k-space as a vector relating to the vector of magnetic field gradient and the signal we acquired in a perticular time is a fourier relationship between the signal in k-space and the spin density, which in this case is a perticular positions in two or three dimensions.

![nmr12]()

Let's look at how we might perform an experiment with a real pulse sequence that enables us to acquire a two-dimensional image.
We start like before with a 90-degree pulse. and now we are going to acquire the signal in reciprocal space in k-space (two-dimension of space), in this case instead of single line of acquisition we are going to have multiple lines represented by a matrix of points (fig nmr13). each of the points in the matrix needs to acquire the signal that will enable us to get a full two-dimensional representation of the signal in the two-dimensions of k-space.
here we are going to carry out a traverse through k-space and look at a way of describing the history of what the spins do with time under the application of magnetic field gradient pulses.
with applying two gradients (one along x-axis another along y-axis) at the same time we can observe the signal for all points in a specific row in the matrix (actually here we use the technic "spin echo" to improve homogeneity).
to fill out all lines of the matrix we must repeat the process for all lines (for each line we use a different areas of pulse along y-axis), this process is known as "phase encoding gradient".
after acquiring the signals they should be transformed using a cartesian fourier transformation.

![nmr13]()


### Acquired Signal and Corresponding Magnetic Resonance Image
![nmr15]()


