package rklt;

import java.util.concurrent.Callable;

import GiciAnalysis.ImageCovariance;
import GiciAnalysis.ImageStatistical;
import GiciException.ParameterException;
import GiciMatrix.MatrixAlgebra;
import GiciTransform.Permutation;
import GiciTransform.ZeroMean;

public class ClusteredRKLT implements Callable {
	
	/* debug/plot function */
	private static void dumpComponentEnergy(float[][][] image) {
		//if (true) return;
		System.out.print("Energy per component: ");
		ImageStatistical is = new ImageStatistical(image);

		double [] d = is.getEnergy(); //is.getVariance(); //

		for (int i = 0; i < d.length; i++) {
			System.out.print(d[i]);

			if (i != d.length - 1) {
				System.out.print(", ");
			} else {
				System.out.print("\n");
			}
		}
	}
	
//	 dump the covariance matrix with a matlab-like format
	public static void dumpImageCovariance(float[][][] image) {
		float[][] cm =  MatrixAlgebra.toComplete(ImageCovariance.generateCovarianceMatrix(image, 0, 0.01f, 0.01f));

		System.out.print("Covariance matrix:\n[");
		for (int i = 0; i < cm.length; i++) {
			System.out.print("[");
			for (int j = 0; j < cm[i].length; j++) {
				System.out.print(cm[i][j] + (j != cm[i].length - 1? ",": ""));
			}
			System.out.print("]" + (i != cm.length - 1? ",": ""));
		}
		System.out.print("]\n");
	}
	
	// Only forward mode is supported, as decoding is universal
	
	final int[] clusters;
	final int clusterMode;
	
	final int size_z;
	final int size_y;
	final int size_x;
	
	final int dimension = 0;
	
	boolean called = false;
	
	final float subsampling;
	
	// Results
	float[] means;
	final ApplyClusters ac;
	
	public ClusteredRKLT(final float[][][] image, final int clusterMode, final int[] clusters,
			final TransformInformation ti, final float subsampling, final int transformFilter) throws ParameterException {
		super();
		
		this.subsampling = subsampling;
		this.ac = new ApplyClusters(image, ti, false, transformFilter);
		
		size_z = image.length;
		size_y = image[0].length;
		size_x = image[0][0].length;

		
		this.clusters = clusters;
		this.clusterMode = clusterMode;
		assert(clusters.length >= 0);
		
		if (clusterMode == 2) {
			assert (clusters.length >= 2);
		} else if (clusterMode == 3) {
			// Check the cluster vector for correctness
			if (clusters.length % 2 != 1) {
				throw new ParameterException("Level 3 clustering requires an odd number of integers as configuration.");
			}
			
			int checkClusterNumber = clusters[0];
			int checkComponentNumber = size_z;
			
			if (checkComponentNumber % checkClusterNumber != 0) {
				throw new ParameterException("The number of components is not divisible by the number of first-level clusters.");
			}
			
			//checkComponentNumber /= checkClusterNumber;
			
			for (int i = 1; i < clusters.length; i += 2) {
				if (checkComponentNumber / checkClusterNumber < clusters[i]) {
					throw new ParameterException("Selected more components than available on the " + ((i + 3) / 2) + "th level.");
				}
				
				checkComponentNumber = checkClusterNumber * clusters[i];
				
				if (checkComponentNumber % clusters[i + 1] != 0) {
					throw new ParameterException("The number of components (" +  checkComponentNumber
							+ ") is not divisible by the number of clusters (" + clusters[i + 1]
							+ ") of the " + ((i + 3) / 2) + "th level .");
				}
				
				checkClusterNumber = clusters[i + 1];
			}
		}
		
	}
	
	public ClusteredRKLT call() throws Exception {
		called = true;
		
		/* This is intented to show the means but not to alter them, as this
		 * is performend within each rklt.
		 */
		ZeroMean zm = new ZeroMean(ac.getImage(), dimension, false);
		//zm.call();
		means = zm.getMeans();

		/* dump initial energy */
		if (false) {
			/* memory hog, two copies of the original image are required */
			zm = new ZeroMean(ac.getImage(), dimension, true);
			float[][][] copy = zm.applyMeans();
			System.out.println("Original");
			dumpComponentEnergy(ac.getImage());
			System.out.println("Zero Means");
			dumpComponentEnergy(copy);
		}

		// Apply the clustering
		if (clusterMode == 0 || clusterMode == 1) {
			int clusterSize = size_z / clusters[0];
			assert(clusterSize * clusters[0] == size_z);
			
			int levels = 1;
			
			if (clusterMode == 1) {
				levels = (int) (Math.log(clusters[0]) / Math.log(2)) + 1;
			}
			
			int levelClusters = clusters[0];
			
			for (int level = 0; level < levels; level++) {

				for (int currentCluster = 0; currentCluster < levelClusters; currentCluster ++) {
					/* permutation */
					Permutation p = new Permutation(size_z);
					int[] permutation = p.getPermutation();
					
					if (currentCluster == 0 && level > 0) {
						/* apply the multilevel permutation */

						for (int i = 0; i < levelClusters * 2; i++) {
							// Move the first half of each cluster to the front
							// clusters are split at different intervals if its size is odd

							int length; // first half length

							if (i % 2 == 0) {
								length = clusterSize / 2;
							} else {
								length = (clusterSize + 1) / 2;
							}

							int src; // where does it cames from
							int dst; // where it goes

							src = i * clusterSize;
							dst = (i / 2) * clusterSize;

							if (i % 2 == 1) {
								// Goes to the second half of that cluster
								dst += clusterSize - length;
							}

							for (int j = 0; j < length; j++) {
								permutation[dst + j] = src + j;
							}

							// Lets go with the second half
							length = clusterSize - length; // If the other half is the longer one this is the short one.

							src = i * clusterSize;
							dst = levelClusters * clusterSize + (i / 2) * clusterSize;

							src += clusterSize - length; // skip the first part

							if (i % 2 == 1) {
								// Goes to the second half of that cluster
								dst += clusterSize - length;
							}

							for (int j = 0; j < length; j++) {
								permutation[dst + j] = src + j;
							}
						}
					}

					p.setPermutation(permutation);
					p.permutationBringTo(clusterSize * currentCluster, 0, clusterSize);
					
					//p.scream();
					
					ac.applyPermutation(p);
					ac.applyCluster(clusterSize, subsampling);
				}

				dumpComponentEnergy(ac.getImage());
				/* reduce */
				levelClusters = levelClusters /2;
			}
		} else if (clusterMode == 2 || clusterMode == 3) {
			// Dynamic/Static cluster size allocation
			int levelNumber = 0;
			int levelClusterCount = clusters[0];			
			int levelClusterSize = size_z / clusters[0];
			assert(levelClusterSize > 0);
			
			int oldLevelClusterCount = -1;
			int oldLevelClusterSize = -1;
			
			int[] pick = null;
			float[] eigenvalues = new float[size_z];
			
			Permutation p = new Permutation(size_z);
			
			while (true) {
				System.out.println("Level: " + levelNumber);
				
				assert (levelClusterCount * levelClusterSize <= size_z);
			
				p.setToIdentity();
				
				if (levelNumber > 0) {
					/* apply the multilevel permutation (use pick) */						
					int acc = 0;
					
					for (int i = 0; i < pick.length; i++) {
						p.permutationBringTo(oldLevelClusterSize * i, acc, pick[i]);
						acc += pick[i];
					}
					
					assert ((acc / levelClusterSize) * levelClusterSize == acc);
				}
				
				ac.applyPermutation(p);
				ac.setSnapshot();
				
				for (int currentCluster = 0; currentCluster < levelClusterCount; currentCluster ++) {
					/* permutation */
					p.setToIdentity();
					p.permutationBringTo(levelClusterSize * currentCluster, 0, levelClusterSize);

					ac.applyPermutation(p);
					float[] r = ac.applyCluster(levelClusterSize, subsampling);
					
					for (int i = 0; i < r.length; i++) {
						eigenvalues[currentCluster * levelClusterSize + i] = r[i];
					}
				}
				
				/* restore previous cluster state.
				 * this is not needed if this is the last cluster, as there is only one.
				 */
				p.setToIdentity();
				for (int i = levelClusterCount - 1; i >= 0 ; i --) {
					p.permutationBringTo(levelClusterSize * i, 0, levelClusterSize);
				}
				
				ac.applyPermutation(p);
				
				/* debug */
				dumpComponentEnergy(ac.getImage());
				
				/* save old number of clusters */
				oldLevelClusterCount = levelClusterCount;
				oldLevelClusterSize = levelClusterSize;
				
				/* determine the new number of clusters */
				if (clusterMode == 2) {
					/* Are we done yet ? */
					if (levelClusterCount <= 1) {
						break;
					}
					
					/* select the components to be further decorrelated */
					pick = ComponentSelection.getComponentsPerCluster(levelClusterSize, levelClusterCount, eigenvalues, ac.getSnapshot(), clusters[1], size_y * size_x);
					
					System.out.print("Components per cluster (cluster size: " + levelClusterSize + ", clusters: " +  levelClusterCount + "): ");
					MatrixAlgebra.printVector(pick);
					
					pick = ComponentSelection.fixComponentsPerCluster(pick, levelClusterSize, levelClusterCount, eigenvalues);
					
					System.out.print("Components per cluster (cluster size: " + levelClusterSize + ") (corrected to %16=0): ");
					MatrixAlgebra.printVector(pick);
					
					levelClusterCount = 0;

					for (int i = 0; i < pick.length; i++) {
						levelClusterCount += pick[i];
					}

					assert ((levelClusterCount / levelClusterSize) * levelClusterSize == levelClusterCount);
					levelClusterCount /= levelClusterSize;
				} else if (clusterMode == 3) {
					System.out.println("Cluster size: " + levelClusterSize + ", clusters: " +  levelClusterCount);
										
					pick = new int[levelClusterCount];
					
					/* dump statistics about selected components */
					if (true) {
						System.out.println("selected components");

						for (int i = 0; i < 8; i++) {					
							pick = ComponentSelection.getComponentsPerCluster(levelClusterSize, levelClusterCount, eigenvalues, ac.getSnapshot(), i, size_y * size_x);
							System.out.print("m: " + i + " - ");
							MatrixAlgebra.printVector(pick);
						}
					}

					/* Are we done yet ? (no more cluster levels in the vector) */
					if (2 * levelNumber + 1 >= clusters.length) {
						break;
					}
					
					/* the number of components to further decorrelated is manually selected */
					levelClusterCount = clusters[2 * (levelNumber + 1)];
									
					for (int i = 0; i < pick.length; i++) {
						pick[i] = clusters[2 * levelNumber + 1];
					}
					
					levelClusterSize = oldLevelClusterCount * clusters[2 * levelNumber + 1] / levelClusterCount;
					
					assert(levelClusterSize > 0);
					
					// check that the number of components selected form the next cluster
					assert(oldLevelClusterCount * clusters[2 * levelNumber + 1] % levelClusterCount == 0);
				}
				
				//assert (levelClusterCount * levelClusterSize == size_z);
				
				/* next level */
				levelNumber++;
			}
		}

		return this;
	}

	/* lossy equivalent getters */
	public float[] getMeans() {
		assert(called);
		return means;
	}

	public float[][] getLossyKLTMatrix() {
		assert(called);
		return ac.getLossyKLTMatrix();
	}
}
