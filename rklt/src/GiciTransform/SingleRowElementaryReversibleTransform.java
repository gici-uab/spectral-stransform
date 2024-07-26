package GiciTransform;

import java.util.concurrent.Callable;

public strictfp class SingleRowElementaryReversibleTransform implements Callable {

	final int[] permutation;
	final float[][] serms;
	long[][][] image;
	boolean forward;
	
	public SingleRowElementaryReversibleTransform(final int[] permutation,
			final float[][] serms, float[][][] image, final boolean forward) {

		this.permutation = permutation;
		this.serms = serms;
		this.forward = forward;
		
		// Asserts
		assert (this.serms.length == image.length + 1);
		
		// copy: no in place calculation
		this.image = new long[image.length][image[0].length][image[0][0].length];
		
		for (int z = 0; z < image.length; z++) {
			for (int y = 0; y < image[0].length; y++) {
				for (int x = 0; x < image[0][0].length; x++) {
					this.image[z][y][x] = (int)image[z][y][x];
				}
			}
		}
		// end of copy
	}

	public SingleRowElementaryReversibleTransform(final SingleRowElementaryReversibleMatrix erm,
			float[][][] image, final boolean forward) {

		this.permutation = erm.getPermutation();
		this.serms = erm.getSingleRowElementarySteps();
		//this.image = image;
		this.forward = forward;
		
		// Asserts 
		assert (this.serms.length == image.length + 1);
		
		// copy: no in place calculation
		this.image = new long[image.length][image[0].length][image[0][0].length];
		
		for (int z = 0; z < image.length; z++) {
			for (int y = 0; y < image[0].length; y++) {
				for (int x = 0; x < image[0][0].length; x++) {
					this.image[z][y][x] = (int)image[z][y][x];
				}
			}
		}
		// end of copy
		
	}

	public float[][][] call() {

		int N = image.length;
		
		//MatrixAlgebra.printMatrix(serms);
		
		if (forward) {
			for (int y = 0; y < image[0].length; y++) {
				for (int x = 0; x < image[0][0].length; x++) {

					for (int k = -1; k < N; k++) {
						int k1 = (k + N) % N; // k=-1 => k=N-1

						double acc = 0;
						
						for (int z = 0; z < N; z++) {
							if (z == k1)
								continue;
							
							//System.out.println("f: " + serms[k+1][z] + " " + image[permutation[z]][y][x] + " " + serms[k+1][z] * image[permutation[z]][y][x]);
							
							acc += (double)serms[k+1][z] * image[z][y][x];
						}

						//System.out.println("r: " + Math.round(acc) + " " + acc);
						assert(!Double.isInfinite(acc) && !	Double.isNaN(acc));
												
						image[k1][y][x] += Math.rint(acc);
						
						// Anything outside this boundary is not lossless
						assert(image[k1][y][x] < ((long)1 << 60) &&
								image[k1][y][x] > -((long)1 << 60));
						
						//System.out.println("ir: " + image[k1][y][x]);
					}
					
					// Apply the permutation
					long[] tl = new long[N];
					
					for (int z = 0; z < N; z++) {
						tl[z] = image[z][y][x];
					}

					for (int z = 0; z < N; z++) {
						image[z][y][x] = tl[permutation[z]];
					}
				}
			}
		} else {
			// Reverse
			
			int[] ipermutation = new int[permutation.length];
			
			for (int i = 0; i < permutation.length; i++) {
				ipermutation[permutation[i]] = i;
			}

			
			for (int y = 0; y < image[0].length; y++) {
				for (int x = 0; x < image[0][0].length; x++) {

					// Undo the permutation
					long[] tl = new long[N];
					
					for (int z = 0; z < N; z++) {
						tl[z] = image[z][y][x];
					}

					for (int z = 0; z < N; z++) {
						image[z][y][x] = tl[ipermutation[z]];
					}
					
					for (int k = N - 1; k >= -1; k--) {
						int k1 = (k + N) % N; // k=-1 => k=N-1

						double acc = 0;
						
						// It is very important here to repeat the operations in
						// the very same order so the floating point unit
						// can deliver the exact same result (although incorrect).
						for (int z = 0; z < N; z++) {
							if (z == k1)
								continue;
							
							//System.out.println("f: " + serms[k+1][z] + " " + image[ipermutation[z]][y][x] + " " + serms[k+1][z] * image[ipermutation[z]][y][x]);
							
							acc += (double)serms[k+1][z] * image[z][y][x];
						}

						//System.out.println("i: " + Math.round(acc) + " " + acc);
						
						image[k1][y][x] -= Math.rint(acc);
						
						//System.out.println("ii: " + image[k1][y][x]);
					}
				}
			}
		}
		
		// copy: no in place calculation
		float[][][] t = new float[image.length][image[0].length][image[0][0].length];
		
		for (int z = 0; z < image.length; z++) {
			for (int y = 0; y < image[0].length; y++) {
				for (int x = 0; x < image[0][0].length; x++) {
					
					// Anything outside this boundary is not lossless
					assert(image[z][y][x] < (1 << 23) &&
							image[z][y][x] > -(1 << 23));
					
					t[z][y][x] = image[z][y][x];
				}
			}
		}
		
		return t;
	}


}
