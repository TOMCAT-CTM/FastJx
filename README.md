# FastJx
FastJx routines for TOMCAT / UKCA / UKESM

This repository contains programs and input data which can be used to average absorption cross sections onto the 18 extended wavelength bins used by the FastJx photolysis scheme, the default scheme in UKCA.

addX_FJX_so3_v2.f  - Fortran program to read input cross sections (in format of JPL table) and write averaged cross sections. This is adapted from original code by Michael Prather

compile - Various example compile lines for above program (Portland, Gnu, Intel)

JPL_data - directory of input cross sections. These are text files in format give by NASA/JPL evaluations.

Out_xsect - Example output of the 'addX_FJX_so3_v2.f' program. These files are cross sections at two temperatures over the 18 wavelength bins and are inserted into input file for the 3D model.
