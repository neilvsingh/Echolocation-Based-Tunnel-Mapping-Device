# Echolocation-Based-Tunnel-Mapping-Device
# Echolocation-Based Tunnel Mapping Device

This project is an echolocation-based tunnel mapping device developed as a coursework project for Fundamentals of Linear Systems. The goal of the system is to use high-frequency acoustic signals to detect nearby surfaces, estimate distances, and generate a two-dimensional map of a tunnel-like environment.

The device concept uses a speaker to transmit short high-frequency pulses and a microphone to receive the reflected echo signals. These return signals are passed through an analog band-pass filter to isolate the useful frequency range and reduce unwanted low-frequency environmental noise. The filtered signal is then sent through an analog-to-digital converter so it can be processed digitally in MATLAB.

The MATLAB program analyzes the received signals using core linear systems techniques, including convolution, Fourier transforms, digital filtering, and phase estimation. Convolution is used to detect the time delay between the transmitted pulse and reflected signal, allowing the system to estimate distance using time-of-flight. Fourier analysis helps examine frequency content and identify changes in the reflected signal that may relate to surface material properties.

The simulation also includes a 2D mapping component that estimates obstacle locations, surface angles, and reflected material characteristics. This project demonstrates how signal processing can be applied to navigation, mapping, and environmental sensing.
