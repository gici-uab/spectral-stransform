/*
 * GICI Library -
 * Copyright (C) 2007  Group on Interactive Coding of Images (GICI)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Group on Interactive Coding of Images (GICI)
 * Department of Information and Communication Engineering
 * Autonomous University of Barcelona
 * 08193 - Bellaterra - Cerdanyola del Valles (Barcelona)
 * Spain
 *
 * http://gici.uab.es
 * gici-info@deic.uab.es
 */
package GiciAnalysis;


/**
 * This class receives an image and calculates some statistical information about the image.
 *
 * @author Group on Interactive Coding of Images (GICI)
 * @version 1.0
 */
public class ImageStatistical{
	
	/**
	 * Array where min (minMax[component][0]) and max (minMax[component][1]) values of each component will be stored.
	 * <p>
	 * All values allowed.
	 */
	double minMax[][] = null;

	/**
	 * Global min (minMax[0]) and max (minMax[1]) values of the image.
	 * <p>
	 * All values allowed.
	 */
	double totalMinMax[] = null;

	/**
	 * Average of each component.
	 * <p>
	 * All values allowed.
	 */
	double average[] = null;

	/**
	 * Average of whole image.
	 * <p>
	 * All values allowed.
	 */
	double totalAverage=0;

	/**
	 * Center of range for each image component.
	 * <p>
	 * All values allowed.
	 */
	double centerRange[] = null;

	/**
	 * Center of range of whole image.
	 * <p>
	 * All values allowed.
	 */
	double totalCenterRange=0;

	/**
	 * Array where are stored the how many times is appeared all values in Byte image.
	 * <p>
	 * All values allowed-
	 */
	int[][] countedValues = null;
	
	/**
	 * Energy of each component.
	 * <p>
	 * All values allowed.
	 */
	double energy[] = null;
	
	/**
	 * Energy of whole image.
	 * <p>
	 * All values allowed.
	 */
	double totalEnergy = 0;
	
	/**
	 * Variance of each component.
	 * <p>
	 * All values allowed.
	 */
	double variance[] = null;
	
	/**
	 * Varianze of whole image.
	 * <p>
	 * All values allowed.
	 */
	double totalVariance = 0;
	
	/**
	 * Histogram of image.
	 * <p>
	 * Only values between the mininum and maximum values of the image are allowed.
	 */
	int hist[][] = null;
	
	/**
	 * Histogram of whole image
	 * <p>
	 * Only values between the mininum and maximum values of the image are allowed.
	 */
	int totalHist[] = null;
	
	/**
	 * Entropy of the image.
	 */
	double entropy[] = null;
	
	/**
	 * Entropy of whole image.
	 */
	double totalEntropy = 0.0D;
	
	/**
	 * Constructor that does all the operations to calculate min and max, average and center range of the image.
	 *
	 * @param imageSamples a 3D float array that contains image samples
	 */
	public ImageStatistical(float[][][] imageSamples){
		//Size set
		int zSize = imageSamples.length;
		int ySize = imageSamples[0].length;
		int xSize = imageSamples[0][0].length;

		//Memory allocation
		minMax = new double[zSize][2];
		average = new double[zSize];
		centerRange = new double[zSize];
		totalMinMax = new double[2];
		countedValues = new int[zSize][256];
		energy = new double[zSize];
		variance = new double[zSize];
		//Initializations
		totalMinMax[0] = Float.POSITIVE_INFINITY;
		totalMinMax[1] = Float.NEGATIVE_INFINITY;
		for(int z = 0; z < zSize; z++){
			minMax[z][0] = Float.POSITIVE_INFINITY;
			minMax[z][1] = Float.NEGATIVE_INFINITY;
			average[z] = 0;
			energy[z] = 0;
			variance[z] = 0;
		}
		
		
		//Calculus
		for(int z = 0; z < zSize; z++){
			for(int y = 0; y < ySize; y++){
				for(int x = 0; x < xSize; x++){
					energy[z] += (imageSamples[z][y][x] * imageSamples[z][y][x]);
					//to count values
					if(Math.abs(imageSamples[z][y][x]) < 256){
						countedValues[z][(int)Math.abs(imageSamples[z][y][x])]++;
					}
					//Min and max
					if(imageSamples[z][y][x] < minMax[z][0]){
						minMax[z][0] = imageSamples[z][y][x];
						if(imageSamples[z][y][x] < totalMinMax[0]){
							totalMinMax[0] = imageSamples[z][y][x];
						}
					}
					if(imageSamples[z][y][x] > minMax[z][1]){
						minMax[z][1] = imageSamples[z][y][x];
						if(imageSamples[z][y][x] > totalMinMax[1]){
							totalMinMax[1] = imageSamples[z][y][x];
						}
					}
					//Average
					average[z] += imageSamples[z][y][x];
				}
			}
			average[z] /= (xSize*ySize);
			
				
			for(int y = 0; y < ySize; y++){
				for(int x = 0; x < xSize; x++){
					variance[z] += ((imageSamples[z][y][x] - average[z]) * (imageSamples[z][y][x] - average[z]));
				}
			}
			variance[z] /= (imageSamples[z].length * imageSamples[z][0].length);
		}
		
		//centerRange, totalVariance
		for(int z = 0; z < zSize; z++){
			centerRange[z] = (minMax[z][0] + minMax[z][1]) / 2;
		}
		
		//Totals
		totalAverage = 0F;
		totalCenterRange = 0F;
		for(int z = 0; z < zSize; z++){
			totalAverage += average[z];
			totalCenterRange += centerRange[z];
			totalEnergy += energy[z];
		}
		totalAverage /= zSize;
		totalCenterRange /= zSize;
		
		// Histogram, totalVariance
		hist = new int[zSize][];
		totalHist = new int[(int)totalMinMax[1] - (int)totalMinMax[0] + 1];
		int totalMin = (int)totalMinMax[0];
		for(int z = 0; z < zSize; z++) {
			int min = (int)minMax[z][0];
			hist[z] = new int[(int)minMax[z][1] - (int)minMax[z][0] + 1];
			for (int i = 0; i < hist[z].length; i++) {
				hist[z][i] = 0;
			}
			for(int y = 0; y < ySize; y++) {
				for(int x = 0; x < xSize; x++){
					totalVariance += (imageSamples[z][y][x] - totalAverage) * (imageSamples[z][y][x] - totalAverage);
					hist[z][(int)imageSamples[z][y][x]-min]++;
					totalHist[(int)imageSamples[z][y][x]-totalMin]++;
				}
			}
		}
		totalVariance /= (imageSamples.length * imageSamples[0].length * imageSamples[0][0].length);
		
		// Entropy
		entropy = new double[zSize];
		for(int z = 0; z < zSize; z++) {
			entropy[z] = 0.0D;
			for (int i = 0; i < hist[z].length; i++) {
				double p = ((double)hist[z][i]) / (ySize*xSize);
				if (p > 0.0D) {
					entropy[z] -= p * Math.log(p);
				}
			}
			entropy[z] /= Math.log(2.0D);		
		}
		
		totalEntropy = 0.0D;
		for (int i = 0; i < totalHist.length; i++) {
			double totalP = ((double)totalHist[i]) / (zSize*ySize*xSize);
			if (totalP > 0.0D) {
				totalEntropy -= totalP * Math.log(totalP);
			}
		}
		totalEntropy /= Math.log(2);

	}

	/**
	 * Constructor that does all the operations to calculate min and max, average and center range of the image.
	 *
	 * @param imageSamples a 3D float array that contains image samples
	 */
	public ImageStatistical(int[][][] imageSamples){
		//Size set
		int zSize = imageSamples.length;
		int ySize = imageSamples[0].length;
		int xSize = imageSamples[0][0].length;

		//Memory allocation
		minMax = new double[zSize][2];
		countedValues = new int[zSize][256];
		average = new double[zSize];
		centerRange = new double[zSize];
		totalMinMax = new double[2];
		double E[] = new double[zSize];
		double totalE = 0;
		//Initializations
		totalMinMax[0] = Float.POSITIVE_INFINITY;
		totalMinMax[1] = Float.NEGATIVE_INFINITY;
		for(int z = 0; z < zSize; z++){
			minMax[z][0] = Float.POSITIVE_INFINITY;
			minMax[z][1] = Float.NEGATIVE_INFINITY;
			average[z] = 0;
			E[z] = 0;
		}
		//Calculus
		for(int z = 0; z < zSize; z++){
			for(int y = 0; y < ySize; y++){
				for(int x = 0; x < xSize; x++){
					energy[z] += (imageSamples[z][y][x] * imageSamples[z][y][x]);
					//to count values
					countedValues[z][(int)imageSamples[z][y][x]]++;
					//Min and max
					if(imageSamples[z][y][x] < minMax[z][0]){
						minMax[z][0] = imageSamples[z][y][x];
						if(imageSamples[z][y][x] < totalMinMax[0]){
							totalMinMax[0] = imageSamples[z][y][x];
						}
					}
					if(imageSamples[z][y][x] > minMax[z][1]){
						minMax[z][1] = imageSamples[z][y][x];
						if(imageSamples[z][y][x] > totalMinMax[1]){
							totalMinMax[1] = imageSamples[z][y][x];
						}
					}
					//Average
					average[z] += (imageSamples[z][y][x] / (xSize*ySize));
				}

			}
			E[z] /= (imageSamples[z].length * imageSamples[z][0].length);
			totalE += E[z];
			
			for(int y = 0; y < ySize; y++){
				for(int x = 0; x < xSize; x++){
					variance[z] = (imageSamples[z][y][x] - E[z]) * (imageSamples[z][y][x] - E[z]);
				}
			}
			variance[z] /= (imageSamples[z].length * imageSamples[z][0].length);
		}
		totalE /= zSize;
		
		//centerRange, totalVariance, hist
		hist = new int[zSize][];
		totalHist = new int[(int)totalMinMax[1] - (int)totalMinMax[0] + 1];
		int totalMin = (int)totalMinMax[0];
		for(int z = 0; z < zSize; z++){
			int min = (int)minMax[z][0];
			hist[z] = new int[(int)minMax[z][1] - (int)minMax[z][0] + 1];
			for (int i = 0; i < hist[z].length; i++) {
				hist[z][i] = 0;
			}
			centerRange[z] = (minMax[z][0] + minMax[z][1]) / 2;
			for(int y = 0; y < ySize; y++){
				for(int x = 0; x < xSize; x++){
					totalVariance = (imageSamples[z][y][x] - totalE) * (imageSamples[z][y][x] - totalE);
					hist[z][imageSamples[z][y][x]-min]++;
					totalHist[(int)imageSamples[z][y][x]-totalMin]++;
				}
			}
		}
		
		// Entropy
		entropy = new double[zSize];
		for(int z = 0; z < zSize; z++) {
			entropy[z] = 0.0D;
			for (int i = 0; i < hist[z].length; i++) {
				double p = ((double)hist[z][i]) / (ySize*xSize);
				if (p > 0.0D) {
					entropy[z] -= p * Math.log(p);
				}
			}
			entropy[z] /= Math.log(2.0D);
		}
		
		totalEntropy = 0.0D;
		for (int i = 0; i < totalHist.length; i++) {
			double totalP = ((double)totalHist[i]) / (zSize*ySize*xSize);
			if (totalP > 0.0D) {
				totalEntropy -= totalP * Math.log(totalP);
			}
		}
		totalEntropy /= Math.log(2);
				
		//Totals
		totalAverage = 0F;
		totalCenterRange = 0F;
		for(int z = 0; z < zSize; z++){
			totalAverage += average[z];
			totalCenterRange += centerRange[z];
			totalEnergy += energy[z];
		}
		totalAverage /= zSize;
		totalCenterRange /= zSize;
		totalVariance /= (imageSamples.length * imageSamples[0].length * imageSamples[0][0].length);
	}
	/**
	 * @return countedValues definition in this class
	 */
	public int[][] getcountedValues(){
		return(countedValues);
	}
	
	/**
	 * @return minMax definition in this class
	 */
	public double[][] getMinMax(){
		return(minMax);
	}

	/**
	 * @return totalMinMax definition in this class
	 */
	public double[] getTotalMinMax(){
		return(totalMinMax);
	}

	/**
	 * @return average definition in this class
	 */
	public double[] getAverage(){
		return(average);
	}

	/**
	 * @return totalAverage definition in this class
	 */
	public double getTotalAverage(){
		return(totalAverage);
	}

	/**
	 * @return centerRange definition in this class
	 */
	public double[] getCenterRange(){
		return(centerRange);
	}

	/**
	 * @return totalCenterRange definition in this class
	 */
	public double getTotalCenterRange(){
		return(totalCenterRange);
	}
	
	/**
	 * @return Energy definition in this class
	 */
	public double[] getEnergy(){
		return(energy);
	}

	/**
	 * @return totalEnergy definition in this class
	 */
	public double getTotalEnergy(){
		return(totalEnergy);
	}
	
	/**
	 * @return Variance definition in this class
	 */
	public double[] getVariance(){
		return(variance);
	}

	/**
	 * @return totalVariance definition in this class
	 */
	public double getTotalVariance(){
		return((float)totalVariance);
	}
	
	/**
	 * @return entropy definition in this class.
	 */
	public double[] getEntropy() {
		return(entropy);
	}
	
	/**
	 * @return totalEntropy definition in this class.
	 */
	public double getTotalEntropy() {
		return(totalEntropy);
	}
}
