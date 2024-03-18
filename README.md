## Data and code for Field Tunable type-I to type-II transition in MoSe2/WS2 heterobilayers

This repository is published alongside **_Electric-field tunable Type-I to Type-II band alignment transition in MoSe2/WS2 heterobilayers_**.\
The data and code contained within can be used to reproduce Figures 2-4 in the main text.

Version 1.0\
Last Updated March 2024

By Jedediah Kistner-Morris\
Quantum Materials Optoelectronics Laboratory\
Department of Physics and Astronomy\
University of California, Riverside, USA

This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.\
If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

###  Data information

The data provided here is used to generate figures 2, 3, and 4 and is separated into PL, dark transport, and photocurrent folders.\
\
**dev1_charge-density-dependent_PL** and **dev1_electric-field-dependent_PL**\
\
_The data in both these folders must be unpacked before it can be used. The PL data is stored in a series of text files. Each file contains a single PL emission spectrum measured at the gate voltage indicated by the file name. Back and Top gates are equal in voltage for the charge dependant data set and oppposite in voltage for the electric field dependant dataset._\
\
**dev2_dark_transport**\
\
_Dark transport data is containted in two files. ...\_pci.dat contains the raw measurement data while ...log.log contains experimental information and measurement metadata._\
\
**dev2_photocurrent**\
\
_Photocurrent data is containted in three files. ...\_pci.npy contains the raw photocurrent data and is a 3d numpy array. ...\_pow.npy contains the laser power data at each measurement point and has the same size and shape of the photocurrent data. ...log.log contains experimental information and measurement metadata._

###  Usage instructions

This code runs on Python (version 3) and requires the following packages:

numpy version 1.23+\
scipy version 1.10+\
matplotlib 3.7+

The python scripts can be run from the command line or from an IDE of your choice.

figure_2.py, figure_3.py, and figure_4.py produce versions of figures 2-4.\
figure_generator.py contains some toolbox functions for generating the figure axes and colormaps.\
extract_interaction_rate_parallel.py extracts the relative rate of interaction by fitting the photocurrent data as described in the main text.\
\
More specific details for each of the scripts are included in the code comments.
