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
# Fig.4) Result