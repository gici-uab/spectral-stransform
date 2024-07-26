package rklt;

public interface ApplyTransform {

	public abstract void apply() throws Exception;
	public abstract void remove() throws Exception;

	public abstract float[][] getPackedData();
	public abstract void setPackedData(final float[][] r);

	public abstract float[][][] getImage();
	public abstract void setImage(float[][][] image);

	public abstract float[][] getKltMatrix();
	public abstract float[] getEigenvalues();

}