package rklt;

import java.util.Arrays;
import java.util.Random;

public class ComponentSelection {

	/* Helper functions for component selection */
	interface itcFunction {
		public float call (float[] eigenvalues, int l, int samples);
	};

	class Energy implements itcFunction {
		public float call (float[] eigenvalues, int l, int samples) {
			double current = 0;
			double total = 0;
			
			for (int s = 0; s < eigenvalues.length; s++) {
				total += eigenvalues[s];
			}

			for (int s = 0; s < l; s++) {
				current += eigenvalues[s];
			}

			final double percent = .99;

			return (float) Math.abs(percent * total - current);
		}
	};
	
	class Akaike implements itcFunction {
		public float call (float[] eigenvalues, int l, int samples) {
			final int m = eigenvalues.length;

			final int remainingValueCount = m - l;

			//float M = l / 2.0f * (remainingValueCount + 1); 
			final float M = l * (2 * m - l);
			
			double acc = 0;
			double avg = 0;

			for (int s = l + 1; s <= m; s++) {
				acc += Math.log(eigenvalues[s - 1]);
				avg += eigenvalues[s - 1];
			}

			avg = avg / (remainingValueCount); 

			double r = -2 * samples * (acc - remainingValueCount * Math.log(avg)) + 2 * M;

			//System.out.println("akaike: " + avg + " " + M + " " + m + "\t" + r);

			return (float) r;
		}
	};

	class MDL implements itcFunction {
		public float call(float[] eigenvalues, int l, int samples) {
			final int m = eigenvalues.length;

			final int remainingValueCount = m - l;

			//float M = l / 2.0f * (remainingValueCount + 1);
			final float M = l / 2.0f * (2 * m - l);
			
			double acc = 0;
			double avg = 0;

			for (int s = l + 1; s <= m; s++) {
				acc += Math.log(eigenvalues[s - 1]);
				avg += eigenvalues[s - 1];
			}

			avg = avg / (remainingValueCount); 

			double r = -2 * samples * (acc - remainingValueCount * Math.log(avg)) + M * Math.log(samples);

			//System.out.println("mdl: " + avg + " " + M + " " + m);
			//System.out.println(r);

			return (float) r;
		}
	};
	
	class EIF implements itcFunction {
		public float call(float[] eigenvalues, int l, int samples) {
			final int m = eigenvalues.length;

			final int remainingValueCount = m - l; 
			
			double avg = 0;

			for (int s = l + 1; s <= m; s++) {
				avg += eigenvalues[s - 1];
			}

			double divisor = (double)remainingValueCount * (double)remainingValueCount * (double)remainingValueCount * (double)samples; 
			avg /= divisor; 

			double r = Math.sqrt(avg);

			return (float) r;
		}
	};
	
	private static int itcMinVectorPos(float[] eigenvalues, itcFunction f, int samples) {
		int minPos = -1;
		float minVal = Float.POSITIVE_INFINITY;
		
		float[] eigenvaluesAbs = Arrays.copyOf(eigenvalues, eigenvalues.length); 
		
		// Same eigenvectors but in some cases in the other direction
		for (int i = 0; i < eigenvalues.length; i++) {
			eigenvaluesAbs[i] = Math.abs(eigenvalues[i]);
		}

		// Quick fix for a correct order 
		for (int i = 0; i < eigenvaluesAbs.length; i++)
			eigenvaluesAbs[i] = -eigenvaluesAbs[i];
		
		Arrays.sort(eigenvaluesAbs);
		
		for (int i = 0; i < eigenvaluesAbs.length; i++)
			eigenvaluesAbs[i] = -eigenvaluesAbs[i];
		
		// Assert descending order
		for (int i = 0; i < eigenvalues.length - 1; i++) {
			assert (eigenvaluesAbs[i] >= eigenvaluesAbs[i+1]);
		}
		
		float[] itc = new float[eigenvaluesAbs.length];
		
		for (int i = 0; i < itc.length; i++) {
			itc[i] = f.call(eigenvaluesAbs, i, samples);
		}
		
		//System.out.print("itc " + f + " ");
		//MatrixAlgebra.printVector(itc);
		
		for (int i = 0; i < itc.length; i++) {
			if (itc[i] < minVal) {
				minPos = i;
				minVal = itc[i];
			}
		}
		
		return minPos;
	}
	
	public static int[] fixComponentsPerCluster(int[] rr, int clusterSize, int clusters, float[] eigenvalues) {
		int[] r = new int[rr.length];
		
		// First, let's fix those clusters that have way to many components
		
		for (int i = 0; i < clusters; i++) {
			r[i] = Math.min(rr[i], clusterSize / 2); 
		}
		
		// we might need some additional components
		
		while (true) {
			int acc = 0;
			
			for (int j = 0; j < r.length; j++)
				acc += r[j];
			
			if (acc % clusterSize == 0)
				break;
			
			// Pick the best next match (biggest next eigenvalue)
			float max = Float.NEGATIVE_INFINITY;
			int maxpos = -1;
			
			for (int i = 0; i < clusters; i++) {
				if (r[i] == clusterSize)
					continue;
				
				if (eigenvalues[i * clusterSize + r[i]] > max) {
					max = eigenvalues[i * clusterSize + r[i]];
					maxpos = i;
				}
			}
			
			//FIXME
			assert (maxpos >= 0);
			r[maxpos]++;
		}
		
		return r;
	}
	
	public static int[] getComponentsPerCluster(int clusterSize, int clusters, float[] eigenvalues, final float[][][] oldImage, final int mode, final int samples) {
		
		//final int size_z;
		final int size_y;
		final int size_x;
		
		/* results is guaranteed to be a multiple of clusterSize */
		int[] r = new int[clusters];

		switch (mode) {
		// Disable unuseful predictors
		case 7:
		case 6:
		case -3:
		case -2:
		case 0: // Half cluster
			for (int i = 0; i < clusters; i++) {
				r[i] = clusterSize / 2;
			}
			
			break;
		case 1: // Average Eigenvalue a.k.a. Latent Root Criterion
			float[] clusterAverage = new float[clusters];
			
			//System.out.print("Eigenvalues: ");
			//MatrixAlgebra.printVector(eigenvalues);
			
			for (int i = 0; i < clusters; i++) {				
				clusterAverage[i] = 0;
				
				int j;
				
				for (j = 0; j < clusterSize; j++) 
					clusterAverage[i] += eigenvalues[i * clusterSize + j];
				
				clusterAverage[i] /= clusterSize;
				
				for (j = 0; j < clusterSize / 2; j++) {
					if (eigenvalues[i * clusterSize + j] < clusterAverage[i])
						break;
				}
				
				r[i] = j;
			}
			
			//System.out.print("Averages: ");
			//MatrixAlgebra.printVector(clusterAverage);
			
			break;
		case 2: // Akaike (AIE)
			for (int i = 0; i < clusters; i++) {				
				float[] eigenValuesCluster = Arrays.copyOfRange(eigenvalues, clusterSize * i, clusterSize * (i + 1));
				r[i] = itcMinVectorPos(eigenValuesCluster, new ComponentSelection().new Akaike(), samples); 
			}
			
			break;
		case 3: // Schwartz and Rissanen (MDL)
			for (int i = 0; i < clusters; i++) {
				float[] eigenValuesCluster = Arrays.copyOfRange(eigenvalues, clusterSize * i, clusterSize * (i + 1));	
				r[i] = itcMinVectorPos(eigenValuesCluster, new ComponentSelection().new MDL(), samples);
			}
			break;
		case 4: // Malinowski (EIF)
			for (int i = 0; i < clusters; i++) {				
				float[] eigenValuesCluster = Arrays.copyOfRange(eigenvalues, clusterSize * i, clusterSize * (i + 1));	
				r[i] = itcMinVectorPos(eigenValuesCluster, new ComponentSelection().new EIF(), samples); 
			}
			
			break;
		case 5: // % of energy/eigenvalue size
			for (int i = 0; i < clusters; i++) {				
				float[] eigenValuesCluster = Arrays.copyOfRange(eigenvalues, clusterSize * i, clusterSize * (i + 1));	
				r[i] = itcMinVectorPos(eigenValuesCluster, new ComponentSelection().new Energy(), samples); 
			}
			break;
		case -6: // HFC
//			for (int i = 0; i < clusters; i++) {				
//				float[][][] oldImageCluster = Arrays.copyOfRange(oldImage, clusterSize * i, clusterSize * (i + 1));	
//				
//				HFC hfc = new HFC(oldImageCluster);
//				
//				r[i] = hfc.call();
//			}
//			
			break;
		case -7: // Parallel Analysis

			//size_z = oldImage.length;
			size_y = oldImage[0].length;
			size_x = oldImage[0][0].length;
			
			for (int i = 0; i < clusters; i++) {				

				float[] eigenValuesCluster = Arrays.copyOfRange(eigenvalues, clusterSize * i, clusterSize * (i + 1));

				// Histogram
//				float[][][] piece = Arrays.copyOfRange(image, i * clusterSize, i * clusterSize + 1); 
//				
//				Reshape rshape = new Reshape(piece, 1, 1, piece[0].length * piece[0][0].length);
//				Histogram hs = (new Histogram(rshape.call()[0][0], 10)).call();
//				
//				MatrixAlgebra.printVector(hs.getBins());
//				MatrixAlgebra.printVector(hs.getCounts());
				
				// Generate some random data
				float[][][] randomData = new float[clusterSize][1][samples];

				Random random = new Random(13);

				//float[][][] oldImage = ac.getSnapshot();
				
				for (int c = 0; c < clusterSize; c++) {
					for (int j = 0; j < samples; j++) {
						int k = random.nextInt(size_y);
						int l = random.nextInt(size_x);

						// This must be before, not after!
						randomData[c][0][j] = oldImage[c + clusterSize * i][k][l];
					}
				}

				ApplyTransform ak = new ApplyRKLT(clusterSize, 0);

				ak.setImage(randomData);

				float[] fakeValues = null;

				try {
					ak.apply();
					fakeValues = ak.getEigenvalues();
				} catch (Exception e) {
					fakeValues = new float[clusterSize];
					Arrays.fill(fakeValues, 0);
					System.out.println("PA FAILED!!!");
				}

				//MatrixAlgebra.printVector(eigenValuesCluster);
				//MatrixAlgebra.printVector(fakeValues);
				
				int c;
				
				for (c = 0; c < clusterSize; c++) {
					if (eigenValuesCluster[c] < fakeValues[c]) {
						break;
					}
				}
				
				r[i] = c;
			}
			break;
			
			default:
				throw new IndexOutOfBoundsException("No such component selection method.");
		}
		
		return r;
	}
}
