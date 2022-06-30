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

from ij import IJ
from ij.gui import Roi, GenericDialog
from ij.io import OpenDialog

from ij.plugin import Stack_Statistics
from ij.plugin.frame import RoiManager

from loci.plugins import BF

import array
import math

#start_time = time.time()

# Part2 - Set up initial parameters
dz = 10                                      #This is the distance between the out of focus image and the focal plane
Lambda = 0.589                               #Wavelength in microns
pi = 3.1415926535
#pi = math.pi
k0 = (2 * pi)/(Lambda)                      #Wavenumber
globaloffset =0.1                            #scaling factor (amplitude) for Hanning bump (alpha in equation (5))
w0 = 2.2                                     #width of the Hanning bump (w0 in equation (5) and (6))
n= 1.34

#Input sample images below. Three images are required, a focal plane image, I0, and two offset images Ia and Ib
#Ia is the ‘above’ image, Ib is the 'below' image, I0 is the in focus image, myBG is the background image 
#Ia = 1.0*(imread('above.tif')) 
#Ib = 1.0*(imread('below.tif')) 
#I0 = 1.0*(imread('focus.tif'))
#BG = 1.0*(imread('bg.tif'))

Ia = IJ.openImage('/Users/guo/Desktop/TIE Code/TIE Code--Python/above.tif')
Ib = IJ.openImage('/Users/guo/Desktop/TIE Code/TIE Code--Python/below.tif')
I0 = IJ.openImage('/Users/guo/Desktop/TIE Code/TIE Code--Python/focus.tif')
BG = IJ.openImage('/Users/guo/Desktop/TIE Code/TIE Code--Python/bg.tif')

Ia_pixels = Ia.getProcessor().getPixels()
Ib_pixels = Ib.getProcessor().getPixels()
I0_pixels = I0.getProcessor().getPixels()
BG_pixels = BG.getProcessor().getPixels()

#Part3 - Transport of intensity equation adapted from Davis
    
#Calculate first derivative of intensity with respect to defocus
#This is equivalent to dI/dz in the paper in equation 1
#Empiracally calculated from measured values so the derivative does not need to be calculated symbolically


dI = array.array('d')
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

wX = array.array('d')
wX.append(-pi + frequencyStepX)
len_wX = 0
while wX[len_wX] < pi + frequencyStepX:
    wX.append(wX[len_wX] + frequencyStepX)
    len_wX += 1

wY = array.array('d')
wY.append(-pi + frequencyStepY)
len_wY = 0
while wY[len_wY] < pi + frequencyStepY:
    wY.append(wY[len_wY] + frequencyStepY)
    len_wY += 1
#print(len_wX)
#print(len_wY)

len_wR = len_wX * len_wY
wR = array.array('d')
for i in range(len_wX):
    for j in range(len_wY):
        wR.append(math.sqrt(wX[i]**2 + wY[j]**2))
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
        hanning[i * len_wX + j] = globaloffset * (1 + math.cos(hanning[i * len_wX + j] * pi / w0))

#Transform dI/dz to the Fourier Domain, use zero padding to 2X dimensions
#This takes dI from line 57 above and obtains the fourier transform, DI, from the three initial images

'''
DI = np.fft.fftshift(np.fft.fft2(dI,[Dimy*Padsize, Dimx*Padsize]))
print("DI[0][0].type:", type(DI[0][0]))
print("DI.shape:", DI.shape)

#Define the Laplacian FD transform pair w/ Hanning bump or global offset
TIEH = (1/ ((4*np.pi**2)*(ww**2+hanning)))               #Hanning bump Figure 2
TIEGO = (1/ ((4*np.pi**2)*(ww**2+globaloffset)))         #Global offset as explained in Figure 2


#compute the auxillary function, psi, as outlined in equation (2)
PSIH = - DI * TIEH                                       #PSIH and PSIGO are the auxillary function in the fourier domain
PSIGO = - DI * TIEGO                            
psiH = np.real(np.fft.ifft2(np.fft.ifftshift(PSIH)))     #Inverse of line 92, transforms the fourier form of the phase
psiGo = np.real(np.fft.ifft2(np.fft.ifftshift(PSIGO)))   #PSIH or PSIGO into the spatial form psiH or psiGo

#To obtain phi, the phase, from psi, the auxillary function we need to divide be scaling factors introduced in the
#fourier transform. Strictly speaking this should be integrated however by taking the difference of the intensity
#we are implicitly integrating, allowing us to equate phi=psi*scaling factors
phiH=psiH/(I0*4*np.pi*k0)                                
phiGo=psiGo/(I0*4*np.pi*k0)
print("I0.shape:", I0.shape)
print("phiH.shape:", phiH.shape)
##################################################################################################

#This next part is just displaying the image and ensuring all values are positive, the TIE calculation is now finished.
#crop the phase image ()
phiH = phiH[0:Dimy,0:Dimx]
phiGo = phiGo[0:Dimy,0:Dimx]

#adjust phase map result so that all values are positive
if (1):
   phiH = phiH - np.min(phiH)
   phiGo = phiGo -np.min(phiGo)





#Part 4 - Refractive index extraction

#Calculate the optical path length from the phase phi, the refractive index here is the media refractive index input in step 1.
OPL=(phiH)*(Lambda/(2*(np.pi)*n))
#Calculate the RI of the image from the optical path length (OPL is obtained in metres and dz is input in um so dz is here multiplied by 10^-6)
RI=OPL/(dz*10**-6)






#Part 5 - Display TIE and RI images, save images, and save csv files.


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
