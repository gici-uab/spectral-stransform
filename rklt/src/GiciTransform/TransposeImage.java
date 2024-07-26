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
package GiciTransform;
import GiciException.*;

/**
 * This class receives an image and transposes it if required.Typically it is used in the following way:<br>
 * constructor(receiveing the image to be transposed)<br>
 * setParameters<br>
 * run<br>
 *  
 * @author Group on Interactive Coding of Images (GICI)
 * @version 1.0
 *
 */
public class TransposeImage{
	
	/**
	 * Original image samples (index meaning [z][y][x]).
	 * <p>
	 * All values allowed.
	 */
	float[][][] imageSamples = null;

	/**
	 * Indicates whether each component should be transposed after reconstruction
	 * <p>
	 * Valid values are:<br>
	 *   <ul>
	 *     <li> 0 - transpose image
	 *     <li> 1 - do not transpose image
	 *   </ul> 
	 */
	int[] transposeImg = null;
	
	/**
	 * To know if parameters are set.
	 * <p>
	 * True indicates that they are set otherwise false.
	 */
	boolean parametersSet;
	
	/**
	 * Constructor of TransposeImage. It receives the image to be transposed.
	 * 
	 * @param imageSamples definition in {@link #imageSamples}
	 */
	public TransposeImage(float[][][] imageSamples){
		//Image data copy
		this.imageSamples = imageSamples;
	}
	
	/**
	 * Set the parameter used to transpose the image and check validity.  If the parameter is not initialized, this class initalizes it with default value.
	 * 
	 * @param transposeImg definition in {@link #transposeImg}
	 * 
	 * @throws ParameterException when an invalid or unrecognized exception is detected it will be thrown
	 */
	public void setParameters(int transposeImg[])throws ParameterException{
		this.transposeImg = transposeImg;
		if (transposeImg!=null){
			if (transposeImg.length!=imageSamples.length){
				throw new ParameterException("parameters are not set properly.");
			}
			for(int z=0;z<transposeImg.length;z++){
				if(transposeImg[z]<0 || transposeImg[z]>1){
					throw new ParameterException("parameters are not set properly.");
				}
			}
		} else {// parameter is not initialized
			this.transposeImg = new int[imageSamples.length];
			for(int z=0;z<transposeImg.length;z++){
				transposeImg[z]=0;
			}
		}
		parametersSet = true;
	}
	
	/**
	 * Transposes the image and returns it.
	 * 
	 * @return a float array containing the transposed image
	 * 
	 * @throws ParameterException when an invalid or unrecognized exception is detected it will be thrown
	 */
	public float[][][] run() throws ParameterException{
		//If parameters are not set run cannot be executed
		if(!parametersSet){
			throw new ParameterException("the image cannot be transposed if parameters are not set.");
		}
		
		int zSize = imageSamples.length;
		for(int z=0;z<zSize;z++){
			if(transposeImg[z]!=0){
				int ySize = imageSamples[z].length;
				int xSize = imageSamples[z][0].length;
				if (ySize!=xSize){// The channel is not squared, then it is necessary to create a new channel
					float transposedChannel[][] = new float[xSize][ySize];
					for(int y=0;y<ySize;y++){
						for(int x=0;x<xSize;x++){
							transposedChannel[x][y] = imageSamples[z][y][x];
						}
					}
					imageSamples[z] = null;
					imageSamples[z] = transposedChannel;
				} else {//The channel is squared
					float temp;
					for(int y=0;y<ySize;y++){
						for(int x=0;x<y;x++){
							temp = imageSamples[z][y][x];
							imageSamples[z][y][x] = imageSamples[z][x][y];
							imageSamples[z][x][y] = temp;
						}
					}
				}
				
			} 
		}
		return imageSamples;
	}
	
	/**
	 * According to the given paramater, determines if the image should be transposed or not. 
	 * 
	 * @param transposeImg definition in {@link #transposeImg}
	 * 
	 * @return a boolean which indicates if the image should be transposed. 
	 */
	public static boolean needTranspose(int[] transposeImg){
		boolean transpose = false;
		if (transposeImg!=null){
			int zSize = transposeImg.length;
			for(int z=0; z<zSize && !transpose; z++){
				if (transposeImg[z]==1){
					transpose = true;
				}
			}
		}
		return transpose;
	}
}
