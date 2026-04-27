# fluorescence-lifetime-sim

## Overview
This project simulates fluorescence lifetime measurements by modeling the convolution of an exponential fluorescence decay with a detector instrument response function (IRF). The goal is to reproduce how real detectors broaden signals and to explore how sampling and detector parameters affect the observed lifetime data.

## Mathematical Model

The measured signal is modeled as the convolution of the fluorescence decay with the detector instrument response function (IRF):

f(t) = h_f(t) * x(t)\
d(t) = h_d(t) * f(t)\
y(t) = g(t)d(t)

where:

- x(t) is the excitation signal
- h_f(t) = exp(-t/τ) is the fluorescence decay
- h_d(t) is the detector IRF (modeled as a Gaussian)
- f(t) is the emitted fluorescence probablility density function (PDF)
- d(t) is the probability density function (PDF) for the detector output
- y(t) is the measured signal

The signal is then sampled at discrete intervals to simulate digitized detector measurements y[n].

## Running the simulation

Clone the repository:

git clone https://github.com/drew-sims/fluorescence-lifetime-sim.git

Open MATLAB and run:

main.m

Please run while within the source directory of the project
