# Coupled-Pendulum

A numerical and multibody simulation of the coupled pendulum system.

## Overview

This project models the dynamics of coupled pendulum (two pendulums coupled by a spring). The system is analyzed using both linear and nonlinear equations of motion. The modal analysis used to characterized the linear system.

## Features

- Equations of motion
- Linearized model
- Modal analysis
- Free-vibration analysis
- Forced-vibration analysis
- Frequency response analysis

## Repository structure

- 'docs/' : Theoretical documentation
- 'matlab/' : MATLAB implementation
- 'simulink/' : Simscape Multibody model

## Requirements

- MATLAB
- Simulink
- Simscape Multibody

## Usage

1. Open 'matlab/CoupledPendulum_main'
2. Modify the system parameters in the Input section.
3. Modify the simulation settings in the Simulation setup section.
4. Run the script.

For the multibody simulation, open 
'simulink/coupled_pendulum.slx'.

## Documentation

See 'docs/coupled_pendulum.pdf' for the derivation and analysis.
