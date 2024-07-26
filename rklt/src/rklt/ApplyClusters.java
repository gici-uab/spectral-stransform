package rklt;

import GiciMatrix.MatrixAlgebra;
import GiciTransform.Permutation;

/**
 * This class will take care that only one permutation is saved before each cluster stage.
 * @author ian
 *
 */
public class ApplyClusters {
	final float[][][] image;
	
	// Snapshots get permutations applied but not clusters
	final float[][][] snapshotImage;
	final boolean enableSnapshots;
	
	final int dimension = 0;
	float[][] lossy;
	float[] means;
	final TransformInformation ti;
	Permutation totalPermutation;
	boolean permutationPending = false;
	
	final int transformFilter;
	
	public ApplyClusters(final float[][][] image, final TransformInformation ti, boolean enableSnapshots, int transformFilter) {
		super();
		this.image = image;
		this.ti = ti;

		assert(image.length > 0 && image[0].length > 0 && image[0][0].length > 0);
		this.lossy = MatrixAlgebra.identityC(image.length);

		this.enableSnapshots = enableSnapshots;
		
		if (enableSnapshots) {
			snapshotImage = new float[image.length][image[0].length][image[0][0].length];
		} else {
			snapshotImage = null;
		}
		
		setSnapshot();
		
		totalPermutation = new Permutation(image.length);
		
		this.transformFilter = transformFilter;
	}
	
	public ApplyClusters(final float[][][] image, final TransformInformation ti) {
		this (image, ti, true, 0);
	}
	
	/* this is needed for pa */
	public void setSnapshot() {
		if (! enableSnapshots) return;
		
		for(int i = 0; i < image.length; i++) {
			for(int j = 0; j < image[0].length; j++) {
				System.arraycopy(image[i][j], 0, snapshotImage[i][j], 0, image[0][0].length);
			}
		}
	}
	
	public float[][][] getSnapshot() {
		return snapshotImage;
	}
	
	/* apply functions */
	public void applyPermutation(final Permutation p) {
		//MatrixAlgebra.printVector(permutation);
		p.permuteInPlace(image, dimension);
		
		if (enableSnapshots) {
			p.permuteInPlace(snapshotImage, dimension);
		}
		
		totalPermutation.appendPermutation(p);
		permutationPending = true;
		
		lossy = MatrixAlgebra.multiplicationCC(p.getPermutationAsMatrix(),lossy);
	}

	public float[] applyCluster(final int clusterSize) throws Exception {
		return applyCluster(clusterSize, 1);
	}

	public float[] applyCluster(final int clusterSize, final float subsampling) throws Exception {
		/* rklt */
		final ApplyTransform transform;
		
		transform = TransformFilterFactory.getTransformFilter(transformFilter, clusterSize, 0, subsampling); 

		float[][][] imageStep = new float[clusterSize][][];
		System.arraycopy(image, 0, imageStep, 0, clusterSize);

		// Memory savings
		for (int i = 0; i < clusterSize; i++) { image[i] = null; }

		transform.setImage(imageStep);
		transform.apply();
		imageStep = transform.getImage();

		System.arraycopy(imageStep, 0, image, 0, clusterSize);

		// Memory savings
		imageStep = null;

		if (ti != null) {
			ti.pushPermutation(totalPermutation.getPermutation());
			ti.pushDataPacket(transform.getPackedData());
			ti.pushTransformType(transformFilter);
		}

		totalPermutation.setToIdentity();
		permutationPending = false;
		
		// Compute the lossy equivalent
		lossy = MatrixAlgebra.multiplicationCC(MatrixAlgebra.padIdentityBR(transform.getKltMatrix(), image.length - clusterSize), lossy);
		
		return transform.getEigenvalues();
	}

	/* Getters */
	
	public float[][][] getImage() {
		return image;
	}

	public float[][] getLossyKLTMatrix() {
		return lossy;
	}
	
	protected void finalize() {
        if (permutationPending && ! totalPermutation.isIdentity()) {
        	System.out.println ("This transform won't be reversible");
        }
	}
}
