import java.awt.image.BufferedImage;
import java.awt.image.PixelGrabber;
import java.io.File;
import java.io.IOException;
import javax.imageio.ImageIO;
import java.lang.Math;

public class TIEJava {

public void run()
{
	// Part2 - Set up initial parameters
	int dz = 10;
	double Lambda = 0.589;
	double pi = 3.1415926535;
	double k0 = (2*pi)/(Lambda);
	double globaloffset =0.1;
	double w0 = 2.2;
	double n= 1.34;
	
	BufferedImage Ia = ImageIO.read(new File("./above.tif"));
	int Ia_width          = Ia.getWidth();
	int Ia_height         = Ia.getHeight();
	
	BufferedImage Ib = ImageIO.read(new File("./below.tif"));
	int Ib_width          = Ib.getWidth();
	int Ib_height         = Ib.getHeight();
	
	BufferedImage I0 = ImageIO.read(new File("./focus.tif"));
	int I0_width          = I0.getWidth();
	int I0_height         = I0.getHeight();
	
	BufferedImage BG = ImageIO.read(new File("./bg.tif"));
	int BG_width          = BG.getWidth();
	int BG_height         = BG.getHeight();
	
	// Part3 - Transport of intensity equation adapted from Davis
	double temp = k0 / (dz * 1.0);
	int temp1 = Ia_width * Ia_height;
	double [] dI = new double [temp1];
	
	PixelGrabber ga = new PixelGrabber(Ia, 0, 0, -1, -1, true);
    ga.grabPixels();
    double [] pixels1 = (double []) ga.getPixels();
    
    PixelGrabber gb = new PixelGrabber(Ib, 0, 0, -1, -1, true);
    gb.grabPixels();
    double [] pixels2 = (double []) gb.getPixels();
	
	for(int i = 0; i < temp1; i++) dI[i] = (pixels1[i] - pixels2[i]) * temp;
	
	int Dimx = Ia_width;
	int Dimy = Ia_height;
	
	int Padsize = 1;
	
	double frequencyStepX = 2 * pi / (Dimx*Padsize);
	double frequencyStepY = 2 * pi / (Dimy*Padsize);
	
	int sizeX = (int)(2 * pi / frequencyStepX);
	int sizeY = (int)(2 * pi / frequencyStepY);
	double [] wX = new double [sizeX];
	double [] wY = new double [sizeY];
	wX[0] = - pi + frequencyStepX;
	for(int i = 0; i < sizeX; i++) wX[i + 1] =  wX[i] + frequencyStepX;
	wY[0] = - pi + frequencyStepY;
	for(int i = 0; i < sizeY; i++) wY[i + 1] =  wY[i] + frequencyStepY;
	
	double [][] wR = new double [sizeX][sizeY];
	for(int i = 0; i < sizeX; i++)
		for(int j = 0; j < sizeY; j++)
			wR[i][j] = Math.sqrt(wX[i] * wX[i] + wY[j] * wY[j]);
	
	double [][] ww = wR.clone();
	
	double temp2 = pi / w0; 
	for(int i = 0; i < sizeX; i++)
		for(int j = 0; j < sizeY; j++)
			if(wR[i][j] <= temp2) wR[i][j] = temp2;
	
	double [][] hanning = wR.clone();
	for(int i = 0; i < sizeX; i++)
		for(int j = 0; j < sizeY; j++)
			hanning[i][j] = globaloffset * (Math.cos(hanning[i][j] * pi / w0) + 1);

	double [][] DI = fftshift(fft2(dI));
	
	double temp3 = 4.0 * pi * pi;
	double[][] TIEH = ww.clone();
	double[][] TIEGO = ww.clone();
	for(int i = 0; i < sizeX; i++)
		for(int j = 0; j < sizeY; j++)
		{
			TIEH[i][j] = 1 / (temp3 * TIEH[i][j] * TIEH[i][j] + hanning[i][j]);
			TIEGO[i][j] = 1 / (temp3 * TIEH[i][j] * TIEH[i][j] + globaloffset);
		}
	
	double [][] PSIH = new double [sizeX][sizeY];
	double [][] PSIGO = new double [sizeX][sizeY];
	for(int i = 0; i < sizeX; i++)
		for(int j = 0; j < sizeY; j++)
		{
			PSIH[i][j] = - DI[i][j] * TIEH[i][j];
			PSIGO[i][j] = - DI[i][j] * TIEGO[i][j];
		}
	
	double [][] psiH = real(ifft2(ifftshift(PSIH)));
	double [][] psiGo = real(ifft2(ifftshift(PSIGO)));
	
	double [][] phiH = psiH.clone();
	double [][] phiGo = psiGo.clone();
	
	if(Boolean.TRUE)
	{
		double min = 99999999;
		for(int i = 0; i < sizeX; i++)
			for(int j = 0; j < sizeY; j++)
				if(min > phiH[i][j]) min = phiH[i][j];
		for(int i = 0; i < sizeX; i++)
			for(int j = 0; j < sizeY; j++) phiH[i][j] -= min;
	
		min = 99999999;
		for(int i = 0; i < sizeX; i++)
			for(int j = 0; j < sizeY; j++)
				if(min > phiGo[i][j]) min = phiGo[i][j];
		for(int i = 0; i < sizeX; i++)
			for(int j = 0; j < sizeY; j++) phiGo[i][j] -= min;
	}
	
	// Part 4 - Refractive index extraction
	double temp4 = Lambda / (2 * pi * n);
	double [][] OPL = phiH.clone();
	for(int i = 0; i < sizeX; i++)
		for(int j = 0; j < sizeY; j++) OPL[i][j] *= temp4;

	double temp5 = dz * 10e-6;
	double [][] RI = OPL.clone();
	for(int i = 0; i < sizeX; i++)
		for(int j = 0; j < sizeY; j++) RI[i][j] /= temp5;
	
}

}
