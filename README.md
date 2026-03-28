# eeg-synchronisation-analysis
Repository containing controlled EEG synchronisation analyses, including MPC dependence on window length, filter bandwidth, frequency mismatch, noise, and application of seizure-related metrics to real data.

This repository contains controlled experiments and analysis code developed as part of my thesis work on phase-based synchronization in epilepsy, with a particular focus on Mean Phase Coherence (MPC).

## Overview

The aim of this repository is to investigate how MPC behaves under controlled conditions and to understand when synchronisation estimates may reflect true coupling versus methodological confounds.

The repository includes experiments on:
- sample size and window length dependence
- filter bandwidth dependence
- frequency mismatch
- noise and signal-to-noise ratio dependence

In addition, it includes MATLAB scripts applying seizure-related metrics to real sqEEG data.

## Relevance

This repository demonstrates practical experience with:
- EEG signal processing
- band-pass filtering
- phase extraction using the Hilbert transform
- sliding-window analysis
- biomarker-oriented analysis in epilepsy
- application of methods to real electrophysiological data

## Repository structure

- `notebooks/` contains controlled experiments investigating MPC behaviour
- `matlab/` contains seizure-related metric scripts applied to real sqEEG recordings
- `figures/` contains example output figures
- `docs/` contains short methodological notes

## Notes

No patient data is included in this repository. Example analyses are based on simulations and code only.

## Author

Madeleine Fairey  
Biomedical Engineering, Universitat Pompeu Fabra
