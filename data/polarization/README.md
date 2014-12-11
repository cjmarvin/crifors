# Polarization files
------------------

Numbers in files are for Nasmyth focus.
Need a conversion factor of 1/1.843.


order | separation, blue end of FSR | central wavelength | separation, red end of FSR 
      | [microns]                   | [microns]          | [microns]
-------------------------------------------------------------------------------------


## Original explanation
Here are the tables you would need to generate the polarization spectrum. Each spectral order splits in two left- and right-polarized versions of itself. The two images of the same order are not parallel by diverge with increasing wavelength. Alexis computed this divergence for all the orders given in the attached files. Each file has four columns: order numbers, separation in microns at the blue end of the free spectral range, in the central wavelength and in the red wavelength. The numbers are for the Nasmyth focus. To convert to detector you would need to divide them by a fixed factor of 1.843. The two images are symmetric relative to the unpolarized order central line but since we are using the 4-hole decker they are shifted from the centered of the slit image up or down (the 2 nodding positions) by 1.25”.
