package rklt;

import java.io.*;
import java.util.Arrays;

import GiciAnalysis.ImageCovariance;
import GiciAnalysis.OrthogonalityTest;
import GiciException.*;
import GiciFile.*;
import GiciMatrix.MatrixAlgebra;
import GiciTransform.Permutation;


/**
 * @author Group on Interactive Coding of Images (GICI)
 * @version 1.0
 */
public class ReversibleKLT {
	/**
	 * Main method of RKLT application. It takes program arguments, loads images and compare them.
	 *
	 * @param args an array of strings that contains program parameters
	 */
	public static void main(String[] args) throws ErrorException{

//		Parse arguments
		ArgsParser parser = null;
		try{
			parser = new ArgsParser(args);
		}catch(Exception e){
			System.err.println("ARGUMENTS ERROR: " +  e.getMessage());
			e.printStackTrace();
			System.out.println("Please report this error (specifying image type and parameters) to: gici-dev@abra.uab.es");
			System.exit(1);
		}

		String inputFile = parser.getInputImageFile();
		String outputFile = parser.getOutputImageFile();
		int[] inputImageGeometry = parser.getInputImageGeometry();
		int[] outputImageGeometry = parser.getOutputImageGeometry();
		int zBegin = parser.getzBegin();

		if (outputFile != null && outputImageGeometry == null) {
			System.out.println("outputImageGeometry is needed is outputFile is specified.");
			System.exit(1);
		}

		
		/*
		 * Some memory debuging
		 */
		byte[] slack = new byte[1024*1024];
		slack[0] = slack[1];
		
		try {
		
		
//			Image load

			if (inputFile == null) {
				System.out.println("At least inputImage is required.");
				System.exit(1);
			}
		
			LoadFile inputImage = null;
			try{
				if(LoadFile.isRaw(inputFile)){				
					//Check parameters of image geometry
					if((inputImageGeometry[0] <= 0) || (inputImageGeometry[1] <= 0) || (inputImageGeometry[2] <= 0)){
						throw new Exception("Image dimensions in \".raw\" or \".img\" data files must be positive (\"-h\" displays help).");
					}
					if((inputImageGeometry[3] < 0) || (inputImageGeometry[3] > 7)){
						throw new Exception("Image type in \".raw\" or \".img\" data must be between 0 to 7 (\"-h\" displays help).");
					}
					if((inputImageGeometry[4] != 0) && (inputImageGeometry[4] != 1)){
						throw new Exception("Image byte order  in \".raw\" or \".img\" data must be 0 or 1 (\"-h\" displays help).");
					}
					if(zBegin == -1){
						inputImage = new LoadFile(inputFile, inputImageGeometry[0], inputImageGeometry[1], inputImageGeometry[2], inputImageGeometry[3], 
								inputImageGeometry[4], false);
					}else{
						inputImage = new LoadFile(inputFile, inputImageGeometry[0], inputImageGeometry[1], inputImageGeometry[2], inputImageGeometry[3], 
								inputImageGeometry[4], false, zBegin);
					}

				}else{
					inputImage = new LoadFile(inputFile);	
				}
			}catch(Exception e){
				System.err.println("IMAGE LOADING ERROR: " + e.getMessage());
				System.exit(2);
			}

			float[][][] image = inputImage.getImage();

			// Dimension
			int dimension = parser.getDimension();

			if(dimension < 0 || dimension > 2 || dimension != 0){
				System.err.println("Not a valid dimension (right now just 0).");
				System.exit(1);
			}

			String transformFile = parser.getTransformFile();
			try {
				if(parser.getDirection() == 1) {
					// Reverse operation
					// Transform Information load
					TransformInformation ti = null;

					FileInputStream fis = new FileInputStream(transformFile);
					ObjectInputStream in = new ObjectInputStream(fis);
					ti = (TransformInformation)in.readObject();
					in.close();

					int steps = ti.getSteps();

					// Undo the klt
					for (int currentCluster = 0; currentCluster < steps; currentCluster ++) {
						/* rklt */
						float[][] packet = ti.popDataPacket();
						int size = packet[0].length;

						ApplyTransform transform = TransformFilterFactory.getTransformFilter(ti.popTransformType(), size, 0, 1);
						transform.setPackedData(packet);

						// This reverse transform can only be applied to images this big
						assert(size <= image.length);

						float[][][] imageStep = Arrays.copyOfRange(image, 0, size); 

						// Memory savings
						for (int i = 0; i < size; i++) { image[i] = null; };

						transform.setImage(imageStep);
						transform.remove();
						imageStep = transform.getImage();

						System.arraycopy(imageStep, 0, image, 0, size);

						// Memory savings					
						imageStep = null;
						transform = null;

						/* permutation */
						int[] permutation = ti.popPermutation();
						Permutation p = new Permutation(permutation);

						p.reversePermutation();
						p.permuteInPlace(image, dimension);
					}

					assert (ti.getSteps() == 0);
				} else {
					// Forward
					// By default one cluster
					int[] clusters = {1};
					int clusterMode = 3;

					if (parser.getEnableClustering()) {
						clusterMode = parser.getClusteringOptions()[0];

						int[] t = parser.getClusteringOptions();
						clusters = Arrays.copyOfRange(t, 1, t.length);
					}

					// Save things in the right place
					TransformInformation ti = new TransformInformation();

					ClusteredRKLT crklt = (new ClusteredRKLT(image, clusterMode, clusters, ti,
							parser.getEnableSubsampling(),parser.getTransformFilter())).call();

					//save transform file
					if (transformFile != null) {
						FileOutputStream fos = new FileOutputStream(transformFile);
						ObjectOutputStream out = new ObjectOutputStream(fos);
						out.writeObject(ti);
						out.close();
					}

					// dump the covariance matrix with a matlab-like format
					if (true) {
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
					
					// save equivalent matrix information
					if (parser.getDumpEquivalentMatrix() != null) {
						System.out.print("Orthogonality test on the equivalent matrix: [");
						float[] r = OrthogonalityTest.test(MatrixAlgebra.transposeC(crklt.getLossyKLTMatrix()));
						System.out.println(r[0] + ", " + r[1] + ", " + r[2] + ", " + r[3] + ", " 
								+ r[4] + ", " + r[5] + ", " + r[6] + ", " + r[7] + "]");
						
						float[][][] i = {MatrixAlgebra.transposeC(crklt.getLossyKLTMatrix())};
						float[][][] j = {{crklt.getMeans()}};
						SaveFile.SaveFileRaw(i, parser.getDumpEquivalentMatrix() + ".klt", 7, 1);
						SaveFile.SaveFileRaw(j, parser.getDumpEquivalentMatrix() + ".means", 7, 1);
					}
				}
			}catch(Exception e){
				e.printStackTrace();
				System.err.println("PROCESSING ERROR: " + e.getMessage());
				System.exit(2);
			}

			//save file
			if (outputFile != null) {

				// Round?
				if (outputImageGeometry[3] < 6) {
					// Maybe here we can implement some nice predictor...

					for (int z = 0; z < image.length; z++) {
						for (int y = 0; y < image[z].length; y++) {
							for (int x = 0; x < image[z][y].length; x++) {
								image[z][y][x] = (float)Math.rint(image[z][y][x]);
							}
						}
					}
				}

				try{
					int[] outputImageGeometry6 = {outputImageGeometry[0], outputImageGeometry[1], outputImageGeometry[2],
							outputImageGeometry[3], outputImageGeometry[4], 0};
					SaveFile.SaveFileByExtension(image,outputFile,outputImageGeometry6);
				} catch(Exception e){
					e.printStackTrace();
					System.err.println("Gici SaveFile ERROR: " + e.getMessage());
					System.exit(4);
				}
			}

		} catch (OutOfMemoryError e) {
			/* recover some memory */
			slack = null;
			System.gc();
			
			/* dump some useful information */
			e.printStackTrace();
			System.err.println("OOM ERROR: " + e.getMessage());
			System.err.println("command line was: " + parser.getCommandLine());
		}

	}
}
