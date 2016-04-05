# Phase Retrival

A Matlab program for phase retrival in atomic holograms, using the method from James' Thesis[^1]

[TOC]

## phaseRetrivalCaller

The main function of the project. It sets some parameters and imports image data.
Its main part is divided into several steps. Each step calls a corresponding function.

## openxlsFigures

This function imports image data from ieee-be-encoded .xls files to 2D matrices.

## RefPhaseRetrival1

This function calculate the phase information of one of the three input images.

## RefPhaseRetrival2

This function calculate the phase of the reference image according to the result of RefPhaseRetrival1.

## SignalPhaseRetrival1

By subtracting the signal and reference image, we got the hologram of the atomic signal.

[^1]: James, Imaging cold atoms with shot-noise and diffraction limited holographic microscopy