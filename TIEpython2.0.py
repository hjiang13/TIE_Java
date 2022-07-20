#Low-cost, open-access refractive index mapping of algal cells using the transport of intensity equation
#Paper available here: (https://www.biorxiv.org/content/10.1101/640755v1)
#Transport of Intesity Equation code
#Adapted from code availabe with 'Transport of Intensity Equation Microscopy for Dynamic Microtubules' by Q. Tyrell Davis
#Davis paper available here:(https://arxiv.org/abs/1707.04139)


#Contents of this code:
#Part1 - Import neccessary libraries
#Part2 - Set up initial parameters
#Part3 - Transport of intensity equation adapted from Davis
#Part4 - Refractive index extraction
#Part5 - Display TIE and RI images, save images, and save csv files.



#Part1 - Import plotting essentials and necessary numpy/scipy #

from ast import Lambda
from fiji.util.gui import GenericDialogPlus
from ij import IJ, ImagePlus
from ij.gui import Roi, GenericDialog
from ij.io import OpenDialog
from ij.process import FloatProcessor
from ij import WindowManager as wm

from ij.plugin import Stack_Statistics, ImageCalculator
from ij.plugin.frame import RoiManager

from loci.plugins import BF

import array
import math

def split_list(alist, wanted_parts=1):
    """Split a list to the given number of parts."""
    length = len(alist)
    # alist[a:b:step] is used to get only a subsection of the list 'alist'.
    # alist[a:b] is the same as [a:b:1].
    # '//' is an integer division.
    # Without 'from __future__ import division' '/' would be an integer division.
    return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts]
             for i in range(wanted_parts) ]

#start_time = time.time()

# Part2 - Set up initial parameters
default_dz = 10                                      #This is the distance between the out of focus image and the focal plane
default_Lambda = 0.589                               #Wavelength in microns
#pi = 3.1415926535
pi = math.pi
k0 = (2 * pi)/(default_Lambda)                      #Wavenumber
default_globaloffset =0.1                            #scaling factor (amplitude) for Hanning bump (alpha in equation (5))
default_w0 = 2.2                                     #width of the Hanning bump (w0 in equation (5) and (6))
default_n= 1.34

#Input sample images below. Three images are required, a focal plane image, I0, and two offset images Ia and Ib
#Ia is the ‘above’ image, Ib is the 'below' image, I0 is the in focus image, myBG is the background image 
#Ia = 1.0*(imread('above.tif')) 
#Ib = 1.0*(imread('below.tif')) 
#I0 = 1.0*(imread('focus.tif'))
#BG = 1.0*(imread('bg.tif'))
guiPlus = GenericDialogPlus("An enhanced GenericDialog")
gui = GenericDialog("TIE analysis")
guiPlus.addNumericField("dz: ", default_dz, 0)
guiPlus.addNumericField("Lambda: ", default_Lambda, 0)
guiPlus.addNumericField("Global offset:", default_globaloffset, 0)
guiPlus.addNumericField("w0:", default_w0, 0)
guiPlus.addNumericField("n: ", default_n, 0)
guiPlus.addFileField("Above Image path:","Above" )
guiPlus.addFileField("Below Image path:","Below" )
guiPlus.addFileField("Focus Image path:","Focus" )
guiPlus.addFileField("BG Image path:","Background" )
gui.showDialog()
guiPlus.showDialog()

if gui.wasOKed():
    dz = gui.getNextNumber()
    Lambda = gui.getNextNumber()
    globaloffset = gui.getNextNumber()
    k0 = (2 * pi)/(Lambda)
    w0 = gui.getNextNumber()
    n = gui.getNextNumber()

    Ia = guiPlus.getNextString()
    Ib = guiPlus.getNextString()
    I0 = guiPlus.getNextString()
    BG = guiPlus.getNextString()

    Ia_pixels = Ia.getProcessor().getPixels()
    Ib_pixels = Ib.getProcessor().getPixels()
    I0_pixels = I0.getProcessor().getPixels()
    BG_pixels = BG.getProcessor().getPixels()

#Part3 - Transport of intensity equation adapted from Davis
    
#Calculate first derivative of intensity with respect to defocus
#This is equivalent to dI/dz in the paper in equation 1
#Empiracally calculated from measured values so the derivative does not need to be calculated symbolically


dI = array.array('f')
#print(len(Ia_pixels))
for i in range(len(Ia_pixels)):
    dI.append((Ia_pixels[i] - Ib_pixels[i]) * k0 / dz)
    
print(dI[0], Ia_pixels[0] - Ib_pixels[0], k0 / dz)


#width = Ia.getProcessor().getWidth()
#print(width)
#dI = (k0  * (Ia_pixels-Ib_pixels)/(dz))

#This section is coordinates and matrix sizes. It converts from the spatial dimensions of the input images into
#frequency units for the Fourier Domain
#Image dimensions

Dimx = Ia.getProcessor().getWidth()
Dimy = Ib.getProcessor().getHeight()
#Dimx = len(dI[0,:])          #dimension of the input images in the x direction
#Dimy = len(dI[:,0])          #dimension of the input images in the y direction
print("Dimx:", Dimx)
print("Dimy:", Dimy)


#set up the coordinates for the Fourier Domain
Padsize = 1                                                                #Increases the size of the image dimensions in the frequency domain
frequencyStepX = 2 * pi / (Dimx*Padsize)                                  #Sets up the size of each step in the Fourier domain X incorporating k
frequencyStepY = 2 * pi / (Dimy*Padsize)                                  #Sets up the size of each step in the Fourier domain Y incorporating k

len_wX = Dimx
wX = array.array('f')
#wX.append(-pi + frequencyStepX)
#len_wX = 0
#while wX[len_wX] < pi + frequencyStepX:
for i in range(1, Dimx + 1):
    wX.append(- pi + i * frequencyStepX)

len_wY = Dimy
wY = array.array('f')
#wY.append(-pi + frequencyStepY)
#while wY[len_wY] < pi + frequencyStepY:
for i in range(1, Dimy + 1):
    wY.append(- pi + i * frequencyStepY)

print("len_wX:", len(wX))
print("len_wY:", len(wY))

len_wR = len_wX * len_wY
wR = array.array('f')
for i in range(len_wY):
    for j in range(len_wX):
        wR.append(math.sqrt(wX[j]**2 + wY[i]**2))
#wX = np.arange(-np.pi+frequencyStepX,np.pi+frequencyStepX,frequencyStepX)  #These variables are fourier domain equivalents to the
#wY = np.arange(-np.pi+frequencyStepY,np.pi+frequencyStepY,frequencyStepY)  #wavenumber k
#WX,WY = np.meshgrid(wX,wY)                                                 #Creates grid with dimensions of wX and Wy
#wR = np.sqrt(WX**2+WY**2)                                                  #wR is kr in equation (4) in the paper            

ww = wR
#print(wR[0])
wR[wR <= pi/w0] = pi/(w0)
#print(wR[0])
#ww = wR                                                     
#wR[wR <= np.pi/w0] = np.pi/(w0)
# wR is used to define the Hanning Bump (set to make hBump value zero outside w0)
# ww is used to calculate the TIE 

################################################################################################

#This is the main section of code, this is where the intensity differential between the three images
#is transformed to the Foruier Domain, where an operation is carried out on it according to the
#Laplacian solution outlined in equation (3) in Davis and where it is transformed back to the spatial domain.

#define the Hanning window, used to stabilize the Fourier Domain TIE
#w0 defines how wide the window is, globaloffset defines the amplitude of the global offset applied

#hanning= globaloffset*(1+np.cos((wR*np.pi)/w0))
hanning = wR
for i in range(len_wX):
    for j in range(len_wY):
        hanning[i * len_wY + j] = globaloffset * (1 + math.cos(hanning[i * len_wY + j] * pi / w0))

#Transform dI/dz to the Fourier Domain, use zero padding to 2X dimensions
#This takes dI from line 57 above and obtains the fourier transform, DI, from the three initial images

#dI.show()
#blank = IJ.createImage("Blank", "32-bit black", Dimy, Dimx, 1)

#blank = IJ.openImage('/Users/guo/Desktop/TIE Code/TIE Code--Python/focus.tif')
#blank_pixels = blank.getProcessor().getPixels()
#blank_pixels1 = array.array("f")
#for i in range(len(blank_pixels)):
#    blank_pixels1.append(float(blank_pixels[i]))
pixel_matrix = split_list(dI, wanted_parts = Dimy)
#print(pixel_matrix)
#pixel_matrix = [list(x) for x in zip(*pixel_matrix)]
DI = ImagePlus("DI", FloatProcessor(pixel_matrix))
#DI.show()
# The input of FloatProcessor is limited to float and int. Short int is invalid.
#IJ.run(DI, "FFT Options...", "complex")
#IJ.run(DI, "FFT", "")
IJ.run(DI, "FFT Options...", "fft complex do");
fft = wm.getImage("FFT of DI")
#fft = wm.getImage("Complex of " + DI.getTitle())
#IJ.run(DI, "FFT 2D 3D", "real=[] imaginary=[] scaling=logarithm");
#fft1 = wm.getImage("Complex of " + DI.getTitle())

#fft_pixels = fft.getProcessor().getPixels()
#print(len(fft_pixels))
#print(fft_pixels)

w_ = fft.getProcessor().getWidth()
h_ = fft.getProcessor().getHeight()
#IJ.makeRectangle(w_/2 - Dimx*Padsize/2, h_/2 - Dimy*Padsize/2, Dimx*Padsize, Dimy*Padsize)
#IJ.run("Duplicate...", " ")
#IJ.makeRectangle(0, 0, Dimy*Padsize, Dimx*Padsize)
#IJ.run("Crop")

#DI = np.fft.fftshift(np.fft.fft2(dI,[Dimy*Padsize, Dimx*Padsize]))
#print("DI[0][0].type:", type(DI[0][0]))
#print("DI.shape:", DI.shape)

#Define the Laplacian FD transform pair w/ Hanning bump or global offset

#TIEH = (1/ ((4*np.pi**2)*(ww**2+hanning)))               #Hanning bump Figure 2
#TIEGO = (1/ ((4*np.pi**2)*(ww**2+globaloffset)))         #Global offset as explained in Figure 2

TIEH = ww
for i in range(len(wR)):
    TIEH[i] = - 1/ ((4 * pi**2)*(TIEH[i]**2 + hanning[i]))

TIEGO = ww
for i in range(len(wR)):
    TIEGO[i] = - 1/ ((4 * pi**2)*(TIEGO[i]**2 + globaloffset))

#compute the auxillary function, psi, as outlined in equation (2)
#PSIH = - DI * TIEH                                       #PSIH and PSIGO are the auxillary function in the fourier domain
#PSIGO = - DI * TIEGO
'''
PSIH = array.array("f")
PSIGO = array.array("f")
for i in range(len(wR)):
    PSIH[i] = - DI[i] * TIEH[i]
for i in range(len(wR)):
    PSIGO[i] = - DI[i] * TIEGO[i]
run("Inverse FFT")                            
psiH = np.real(np.fft.ifft2(np.fft.ifftshift(PSIH)))     #Inverse of line 92, transforms the fourier form of the phase
psiGo = np.real(np.fft.ifft2(np.fft.ifftshift(PSIGO)))   #PSIH or PSIGO into the spatial form psiH or psiGo
'''

Dim_padded = 256
TIEH1 = array.array("f")
TIEGO1 = array.array("f")
for i in range(0, Dim_padded**2):
    TIEH1.append(0.0)
    TIEGO1.append(0.0)
    
for i in range(h_/2 - Dimy*Padsize/2, h_/2 + Dimy*Padsize/2):
    for j in range(w_/2 - Dimx*Padsize/2, w_/2 + Dimx*Padsize/2):
        TIEH1[i * Dim_padded + j] = TIEH[(i - (h_/2 - Dimy*Padsize/2)) * Dimx*Padsize + (j - w_/2 - Dimx*Padsize/2)]
        TIEGO1[i * Dim_padded + j] = TIEGO[(i - (h_/2 - Dimy*Padsize/2)) * Dimx*Padsize + (j - w_/2 - Dimx*Padsize/2)]

pixel_matrix1 = split_list(TIEH1, wanted_parts = Dim_padded)
TIEH_IM = ImagePlus("TIEH", FloatProcessor(pixel_matrix1))
pixel_matrix2 = split_list(TIEGO1, wanted_parts = Dim_padded)
TIEGO_IM = ImagePlus("TIEGO", FloatProcessor(pixel_matrix2))

TIEH_IM.show()
TIEGO_IM.show()

ic = ImageCalculator()
#imp1 = wm.getImage("Complex of DI-1") # Dumplicated Image
imp1 = wm.getImage("Complex of DI")
imp1.show()
imp2 = wm.getImage("TIEH")
#imp2.show()
imp3 = wm.getImage("TIEGO")
res1 = ic.run("Multiply create 32-bit stack", imp1, imp2)
res2 = ic.run("Multiply create 32-bit stack", imp1, imp3)
res1.show()
res2.show()
#imp1 = ic.run("Multiply create 32-bit stack", "Complex of DI", "TIEH")
#imp2 = ic.run("Multiply create 32-bit stack", "Complex of DI", "TIEGO")
#imp1.show()
#imp2.show()
#imp1 = ic.run("Multiply create 32-bit", fft, TIEH_IM)
#imp2 = ic.run("Multiply create 32-bit", fft, TIEGO_IM)

#IJ.run("FD Math...", "image1=[FFT of DI] operation=Convolve image2=[TIEH] result=Result1")

print("1")

IJ.run(fft, "Inverse FFT", "")

IJ.run(res1, "Inverse FFT", "")
IJ.setSlice(2)
IJ.run("Delete Slice")
IJ.makeRectangle(0, 0, Dimy*Padsize, Dimx*Padsize)
IJ.run("Duplicate...", " ")
IJ.run("Rotate 90 Degrees Left")
IJ.run("Flip Vertically")

IJ.run(res2, "Inverse FFT", "")
IJ.setSlice(2)
IJ.run("Delete Slice")
IJ.makeRectangle(0, 0, Dimy*Padsize, Dimx*Padsize)
IJ.run("Duplicate...", " ")
IJ.run("Rotate 90 Degrees Left")
IJ.run("Flip Vertically")

print("2")

imp_filtered1 = wm.getImage("Complex of DI-1-1")
imp_filtered2 = wm.getImage("Complex of DI-2-1")
'''
IJ.setSlice(2)
IJ.run("Delete Slice")
IJ.makeRectangle(0, 0, Dimy*Padsize, Dimx*Padsize)
IJ.run("Crop")
IJ.run("Flip Horizontally")
IJ.run("Flip Vertically")
'''

#imp_pixels = imp_filtered1.getProcessor().getPixels()

psiH = imp_filtered1.getProcessor().getPixels()
psiGo = imp_filtered2.getProcessor().getPixels()

#To obtain phi, the phase, from psi, the auxillary function we need to divide be scaling factors introduced in the
#fourier transform. Strictly speaking this should be integrated however by taking the difference of the intensity
#we are implicitly integrating, allowing us to equate phi=psi*scaling factors
temp = 4 * pi * k0
phiH = array.array('f')
phiGo = array.array('f')

print("len I0_pixels:", len(I0_pixels))
print("len phiH", len(psiH))

for i in range(len(psiH)):
    phiH.append(psiH[i])
    phiGo.append(psiGo[i])

for i in range(len(psiH)):
    phiH[i] = phiH[i] / (I0_pixels[i] * temp)

for i in range(len(psiGo)):
    phiGo[i] = phiGo[i] / (I0_pixels[i] * temp) 

#phiH=psiH/(I0*4*np.pi*k0)                                
#phiGo=psiGo/(I0*4*np.pi*k0)
#print("I0.shape:", I0.shape)
#print("phiH.shape:", phiH.shape)
##################################################################################################

#This next part is just displaying the image and ensuring all values are positive, the TIE calculation is now finished.
#crop the phase image ()

#phiH = phiH[0:Dimy,0:Dimx]
#phiGo = phiGo[0:Dimy,0:Dimx]


#adjust phase map result so that all values are positive
if (1):
#   phiH = phiH - np.min(phiH)
#   phiGo = phiGo -np.min(phiGo)
    min = 99999
    for i in range(len(phiH)):
        if phiH[i] < min:
            min = phiH[i]
    for i in range(len(phiH)):
        phiH[i] -= min

    min = 99999
    for i in range(len(phiGo)):
        if phiGo[i] < min:
            min = phiGo[i]
    for i in range(len(phiGo)):
        phiGo[i] -= min




#Part 4 - Refractive index extraction

#Calculate the optical path length from the phase phi, the refractive index here is the media refractive index input in step 1.
temp = Lambda/(2 * pi * n)
OPL = phiH
for i in range(len(OPL)):
    OPL[i] *= temp
#OPL=(phiH)*(Lambda/(2*(np.pi)*n))
#Calculate the RI of the image from the optical path length (OPL is obtained in metres and dz is input in um so dz is here multiplied by 10^-6)
#RI=OPL/(dz*10**-6)
temp1 = dz*10**-6
RI = OPL
for i in range(len(RI)):
    RI[i] /= temp1




#Part 5 - Display TIE and RI images, save images, and save csv files.

pixel_matrix01 = split_list(OPL, wanted_parts = Dimy)
TIE = ImagePlus("TIE", FloatProcessor(pixel_matrix01))
TIE.show()
IJ.run("Rotate 90 Degrees Left")
IJ.run("Flip Vertically")

pixel_matrix02 = split_list(OPL, wanted_parts = Dimy)
RefractiveIndex = ImagePlus("Refractive Index", FloatProcessor(pixel_matrix02))
RefractiveIndex.show()
IJ.run("Rotate 90 Degrees Left")
IJ.run("Flip Vertically")

'''
plt.figure(figsize=(5,5))           #Change figure size using
plt.title('TIE')                    #Change figure title
plt.imshow(phiH, cmap='gray')       #Specify variable to plot and colormap using
plt.axis('off')                     #Show or hide axes using
#Repeat for other variables, in this case Refractive index
plt.figure(figsize=(5,5)),plt.title('Refractive Index'),plt.imshow(RI, cmap='gray'),plt.axis('off')
plt.show()                          #Show all images


#save images
im = Image.fromarray(phiH);im.save('TIE.tif')
im = Image.fromarray(RI);im.save('RI.tif')

#save data as csv
np.savetxt('TIE.csv', phiH, fmt='%10.8f', delimiter=',')
np.savetxt('RI.csv', RI, fmt='%10.8f', delimiter=',')


print("%s seconds" % (time.time() - start_time))
'''