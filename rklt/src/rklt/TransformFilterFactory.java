package rklt;

public class TransformFilterFactory {
	
	public static ApplyTransform getTransformFilter(final int filter, final int clusterSize,
			final int dimension, final float subsampling) {
		
		ApplyTransform r = null;
		
		switch (filter) {
		case 0:
			r = new ApplyRKLT(clusterSize, 0, subsampling);
			break;
		case 1:
			r = new ApplyKLT(clusterSize, 0, subsampling);
			break;
		default:
			throw new ArrayIndexOutOfBoundsException("Invalid transform");
		}
		
		return r;
	}
}

	
