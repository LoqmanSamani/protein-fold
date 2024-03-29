## Fourier Transformation

In physics, engineering and mathematics, the Fourier transform (FT) is an integral transform that converts a function into a form that describes the frequencies present in the original function. The output of the transform is a complex-valued function of frequency. The term Fourier transform refers to both this complex-valued function and the mathematical operation. When a distinction needs to be made the Fourier transform is sometimes called the frequency domain representation of the original function. 



### Mathematical formulation of the fourier transform


------------------------------------------------------------------------------------------------------------------------------------------


    Fig.1) The fourier transform is a mathematical process which turns a time-domain signal (FID) into a frequency-domain signal (spectrum).


    Intensity at f HZ = area under [FID * cos(2 * pi * f * t)]

    Intensity at f HZ = S_spectrum (f)
    S_FID(t) = FID
    the area under a function is the same as its integral

    
     S_spectrum(f) = ∫ S_FID(t) * cos(2 * pi * f * t) dt 

        integral from 0 to +inf however in practice, as the FID eventually
        decays to zero, we only need to compute the integral from zero up 
        to some time at which the FID is negligible.

        

        On the spectrometer the FID is represented by a set of points which
        have been sampled at regular intervals. To Fourier transform such data,
        rather than computing the integral we multiply each point in the FID by the
        value of the trial cosine wave computed for the time corresponding to that
        point. This gives the product function as a series of data points, and these
        are then summed to give the intensity in the spectrum at the frequency of
        the trial cosine wave.


    S_spectrom(f) = sum(S_FID(ti) * cos(2 * pi * f * ti))

         ti = the time corresponding to the ith data point
         S F1D(ti) = the value of the FID at this time
         N = the total number of data points





    Representing the FID:

    the x- and y—components of the magnetization generated by a 90°(x) pulse can be written:

    1 . My = -M0 * cos(omega * t),     Mx = M0 *  sin(omega * t)

          M0 = the equilibrium magnetization
          omega = the offset (in rad s⁻1)


          If our pulse—acquire experiment had used a 90° pulse about y, rather
          than about x, the equilibrium magnetization would have been rotated onto
          the x-axis:

    2 . Mx = M0 * cos(omega * t),    My = M0 *  sin(omega * t)


         The spectrometer is capable of detecting simultaneously both the
         x- and y—components of the magnetization, each giving rise to 
         separate signals which we will denote Sx and Sy.


    3 . Sx = S0 * cos(omega * t),    Sy = S0 *  sin(omega * t)


          Finally, we need to recognize that the magnetization. and hence the signal.
          will decay over time. We model this by assuming that the signal decays
          exponentially:

    4 .  Sx = S0 * cos(omega * t) * exp(-t/T2),    Sy = S0 *  sin(omega * t) * exp(-t/T2)

          T2 = time constant which characterizes the decay




          Rather than dealing with the x- and y—components separately, it is 
          convenient to bring them together as a complex signal, with the 
          x-component becoming the real part and the y-component the imaginary part.


    5 .  S(t) = Sx + i*Sy

              = S0 * cos(omega * t) * exp(-t/T2) + i[S0 *  sin(omega * t) * exp(-t/T2)]

              = S0 * (cos(omega * t) + sin(omega * t)) * exp(-t/T2)
         
                (cos(teta) + i sin(teta) = exp(i * teta))

              = S0 * exp(i * omega * t) * exp(-t/T2)

                if R = 1/T2   (R in units of s-¹ or Hz)

              = S0 * exp(i * omega * t) * exp(-R * t)

    
              
              If there are several resonances present, then the complex
              time-domain signal is a sum of terms:


    6 . S0 * exp(i * omega * t) * exp(-t/T2) + S0 * exp(i * omega * t) * exp(-t/T2) + ...

         where each resonance i has its own intensity(S0), frequency(omega), and decay constant(T2)


    
         Fourier transformation gives a complex frequency-domain signal or spectrum.
         Normally, the software on the spectrometer only displays the real part of this
         complex spectrum, but it is important to realize that the imaginary part exists,
         even if it is not displayed.



    Fig.2) Fourier transformation of an exponentially decaying time-domain signal.


    7 .  
         S(t) --FT--> S(ω)

         S0 * exp(iΩt) * exp(-Rt) --FT--> [S0 * R / R² + (ω - Ω)²] + i(-S0(ω - Ω)/R² + (ω - Ω))
                                          --------real------------   -------imaginary----------


         The factor of S0 is just an overall scaling, so it is usual to
         deﬁne the Lorentzian absorption and dispersion mode lineshape functions,
         A(ω) and D(ω). without this factor.


    8 . A(ω) = R / R² + (ω - Ω)²,    D(ω) =  i(-(ω - Ω)/R² + (ω - Ω))


        Using the deﬁnitions of the lineshape functions A(ω) and D(ω)

    9 . 
         S(t) --FT--> S(ω)

         S0 * exp(iΩt) * exp(-Rt) --FT--> S0 [A(ω) + i * D(ω)]

    
         As a result, the time-domain signal measured by the spectrometer has an
         essentially arbitrary (and usually unknown) phase Ø associated with it.
     

    10 . S(t) = S0 * exp(iΩt) * exp(-Rt) * exp(i * Ø)

         S0 * exp(iΩt) * exp(-Rt) * exp(i * Ø)   --FT-->   S0 [A(ω) + i * D(ω)] * exp(i * Ø)

   
    11 .  S0[cos(Ø) * A(ω) - sin(Ø) * A(ω)]  +  i*S0[cos(Ø) * D(ω) + sin(Ø) * A(ω)]
          --------------real---------------     ------------imaginary--------------





    Fig.3) Depiction of the effect of a phase shift on the spectrum.

------------------------------------------------------------------------------------------------------------------------


### Fig.1

--------
![fig.1](https://github.com/LoqmanSamani/protein_structure_analysis/blob/systembiology/images/ft1.png)
---------

### Fig.2

----------
![fig.2](https://github.com/LoqmanSamani/protein_structure_analysis/blob/systembiology/images/ft2.png)
----------

### Fig.3

----------
![fig.3](https://github.com/LoqmanSamani/protein_structure_analysis/blob/systembiology/images/ft3.png)
----------


### Implementation of Fourier Transformation in python

```python

from machine_learning.linear_algebra import intro_numpy as np
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
# Fig.4) Result
```

### Fig.4

----------
![fig.4](https://github.com/LoqmanSamani/protein_structure_analysis/blob/systembiology/images/ft4.png)
----------



Reference: ![understanding nmr spectroscopy](https://books.google.de/books?hl=en&lr=&id=PKQlfaK4COoC&oi=fnd&pg=PR17&dq=understanding+nmr+spectroscopy+by+james+keeler&ots=ycPrHp5V8K&sig=zpf_lmOM2uwH_OYzCvYBeMliuLw&redir_esc=y#v=onepage&q=understanding%20nmr%20spectroscopy%20by%20james%20keeler&f=false) by james keeler


Loghman Samani 

November 10, 2023