package rklt;

/**
 * Arguments parser for ffc. This class analyses a string of arguments and extract and check its validity.
 * Usage example:<br>
 * &nbsp; construct<br>
 * &nbsp; [showArgsInfo]<br>
 * &nbsp; [get functions]<br>
 *
 * @author Group on Interactive Coding of Images (GICI)
 * @version 1.0
 */
public class ArgsParser{

	/**
	 * Arguments specificiation. The array describes argument, explain what is used and its default parameters. First index of array is argument; second specifies:<br>
	 *   <ul>
	 *     <li> 0 - short argument specification (i.e. "-i")
	 *     <li> 1 - long argument specification (i.e. "--inputFile")
	 *     <li> 2 - parsing specification of argument ({} indicates mandatority, [] optionality)
	 *     <li> 3 - default values
	 *     <li> 4 - mandatory argument ("1") or non mandatory argument ("0")
	 *     <li> 5 - explanation
	 *   </ul>
	 * <p>
	 * String arguments.
	 */
	String[][] argsSpecification = {
		{"-h", "--help", "", "", "0",
			"Displays this help and exits program."
		},
		{"-i", "--inputImageFile", "{string}", "", "0",
			"Input image. Valid formats are: pgm, ppm, pbm, jpg, tiff, png, bmp, gif, fpx. If image is raw data file extension must be \".raw\" or \".img\" and \"-ig\" parameter is mandatory."
		},
		{"-ig", "--inputImageGeometry", "{int int int int boolean}", "", "0",
			"Geometry of input raw image data. Parameters are:\n    1- zSize (number of image components)\n    2- ySize (image height)\n    3- xSize (image width)\n    4- data type. Possible values are:\n \t 0- boolean (1 byte)\n \t 1- unsigned int (1 byte)\n \t 2- unsigned int (2 bytes)\n \t 3- signed int (2 bytes)\n \t 4- signed int (4 bytes)\n \t 5- signed int (8 bytes)\n \t 6- float (4 bytes)\n \t 7- double (8 bytes)\n    5- Byte order (0 if BIG ENDIAN, 1 if LITTLE ENDIAN)."
		},
		{"-o", "--outputImageFile", "{string}", "", "0",
			"Output image file name. Valid formats are: pgm, ppm, pbm, jpg, tiff, png, bmp, gif, fpx. If image is raw data file extension must be \".raw\" or \".img\" and \"-og\" parameter is mandatory."
		},
		{"-og", "--outputImageGeometry", "{int int int int boolean}", "", "0",
			"Geometry of output raw image data. Parameters are:\n    1- zSize (number of image components)\n    2- ySize (image height)\n    3- xSize (image width)\n    4- data type. Possible values are:\n \t 0- boolean (1 byte)\n \t 1- unsigned int (1 byte)\n \t 2- unsigned int (2 bytes)\n \t 3- signed int (2 bytes)\n \t 4- signed int (4 bytes)\n \t 5- signed int (8 bytes)\n \t 6- float (4 bytes)\n \t 7- double (8 bytes)\n    5- Byte order (0 if BIG ENDIAN, 1 if LITTLE ENDIAN)."
		},
		{"-ca", "--coefficientsApproximation", "{int[ int[ int[ ...]]]}", "0", "0",
			"This parameter specifies the approximation to be applied to coefficients of each channel. First value is for the first channel, second value for the second channel and so on. If only one value is specified, it will be the same for all channels. Valid values are:\n 0 - None  \n 1 - Round  \n 2 - Floor  \n 3 - Ceil "	
		},
		{"-mn", "--minValue", "{float[ float[ float[ ...]]]}", "null", "0",
			"This parameter specifies the mininum value that can appear in a channel. First value is for the first channel, second value for the second channel and so on. If only one value is specified, it will be the same for all channels. If no value is given, no minimum value is considered."	
		},
		{"-mx", "--maxValue", "{float[ float[ float[ ...]]]}", "null", "0",
			"This parameter specifies the maximum value that can appear in a channel. First value is for the first channel, second value for the second channel and so on. If only one value is specified, it will be the same for all channels. If no value is given, no maximum value is considered."	
		},
		{"-ti", "--TransformFile", "{string}", "", "0",
			"Transform file. Usefull to reverse the transform!"
		},
		{"-D", "--dimension", "{int}", "0", "0",
			"Dimension on which the transform will be applied (z=0, y=1, x=2)."
		},
		{"-d", "--direction", "{int}", "", "1",
			"0 for Forward, 1 for Reverse (matrixFile is writen in forward mode and read in reverse)."
		},
		{"-zb", "--zBegin", "{int}", "", "0",
			"First component where the raw file is begun to be load."
		},
		{"-ec", "--enableClustering", "{int int[ int[ int[ ...]]]}", "", "0",
			"Enables clustering (first parameter: mode. second. mode parameter). " +
			"Mode 0 is simple clustering and mode 1 is multilevel clustering. " +
			"Mode 0 and 1 take one extra parameter that is the number of clusters at the first level of clustering. " +
			"Modes 2 and 3 correspond to static and dynamic as defined in the second reference. " + 
			"Mode 2 takes a parameter list: [(#c_1,Th(c_1)),(#c_2,Th(c_2)),...,(#c_L,Th(c_l))], as defined in the second reference. " +
			"Mode 3 takes two parameters: the number of clusters at the first level of clustering, and the eigen-thresholding method to use (1=AE, 4=EIF). " +
			"References: \"Ian Blanes and Joan Serra-Sagristà. Clustered Reversible-KLT for Progressive Lossy-to-Lossless 3d Image Coding. In Data Compression Conference 2009. IEEE Press, 2009\". \"Ian Blanes and Joan Serra-Sagristà. Cost and scalability improvements to the Karhunen-Loêve transform for remote-sensing image coding. IEEE Trans. Geosci. Remote Sens., 2010\"."
		},
		{"-dm", "--dumpEquivalentMatrix", "{string}", "", "0",
			"Matrix transform file. A square matrix for the equivalent lossy transform. (doubles, dimension*dimension)"
		},
		{"-es", "--enableSubsampling", "{float}", "", "0",
			"Enables covariance matrix subsampling (0 > factor >= 1). factor == 1 means no subsampling. " +
			"Factors > 0.1 are not recomended as random selection might sample the same pixel more than once."
		},
		// Kept for compatibility
		{"-lm", "--lossy", "", "", "0",
			"Enables lossy calculations of the KLT."
		},
		{"-t", "--transform", "", "0", "0",
			"Selects the transform to be used inside clusters.\n\t0 - RKLT (default)\n\t1 - KLT (same as --lossy)"
		},
	};

	//ARGUMENTS VARIABLES
	String inputImageFile = null;
	String transformFile = null;
	String outputImageFile = null;
	int[] inputImageGeometry = null;
	int[] outputImageGeometry = null;
	int[] coefficientsApproximation = null;
	float[] minValue = null;
	float[] maxValue = null;
	int dimension = 0;
	int direction = 0;
	int zBegin = -1;
	boolean enableClustering = false; 
	String dumpEquivalentMatrix = null;
	int[] clusteringOptions = null;
	float subsampling = 1;
	int[] trackDependencies = null;
	String trackMaskInput  = null;
	String trackMaskOutput  = null;
	int[] trackMaskGeometry  = null;
	int transform = 0;
	
	String commandLine;
	
	 /**
	  * Class constructor that receives the arguments string and initializes all the arguments
	  * 
	  * @param args the array of strings passed at the command line
	  * 
	  * @throws Exception when an invalid parsing is detected or some problem with method invocation occurs
	  */
	public ArgsParser(String[] args) throws Exception{
		for (String a: args) {
			commandLine += " " + a;
		}
		
		int argNum = 0;
		boolean[] argsFound = new boolean[argsSpecification.length];

		//Arguments parsing
		for(int i = 0; i < argsSpecification.length; i++){
			argsFound[i] = false;
		}
		while(argNum < args.length){
			int argFound = argFind(args[argNum]);
			if(argFound != -1){
				if(!argsFound[argFound]){
					argsFound[argFound] = true;
					int argOptions = argNum + 1;
					while(argOptions < args.length){
						if(argFind(args[argOptions]) != -1){
							break;
						}else{
							argOptions++;
						}
					}
					int numOptions = argOptions - argNum;
					String[] options = new String[numOptions];
					System.arraycopy(args, argNum, options, 0, numOptions);
					argNum = argOptions;
					switch(argFound){
					case 0: //-h  --help
						showArgsInfo();
						System.exit(1);
						break;
					case 1: //-i  --inputImageFile
						inputImageFile = parseString(options);
						if(inputImageFile.endsWith(".raw")){
							argsSpecification[2][4] = "1";
						}
						break;
					case 2: //-ig  --inputImageGeometry
						inputImageGeometry = parseIntegerArray(options, 5);
						break;
					case 3: //-o  --outputImageFile
						outputImageFile = parseString(options);
						break;
					case 4: //-og  --outputImageGeometry
						outputImageGeometry = parseIntegerArray(options, 5);
						break;	
					case 5: //-ca --coefficientsApproximation
						coefficientsApproximation =	parseIntegerArray(options);
						break;
					case 6: //"-mn", "--minValue"
						minValue = parseFloatArray(options);
						break;
					case 7: //"-mx", "--maxValue"
						maxValue = parseFloatArray(options);
						break;
					case 8: //"-ti", "--TransformFile"
						transformFile = parseString(options);
						break;
					case 9: //"-D", "--dimension"
						dimension = parseIntegerPositive(options);
						break;
					case 10: //"-D", "--direction"
						direction = parseIntegerPositive(options);
						break;
					case 11: //"-zb", "--zBegin"
						zBegin = parseIntegerPositive(options);
						break;
					case 12: //"-ec", "--enableClustering"
						enableClustering = true;
						clusteringOptions = parseIntegerArray(options);
						break;
					case 13: //"-dm", "--dumpEquivalentMatrix"
						dumpEquivalentMatrix = parseString(options);
						break;
					case 14: //"-es", "--enableSubsampling"
						subsampling = parseFloatArray(options, 1)[0];
						break;
					case 15: //"-lm", "--lossy"
						transform = 1;
						break;
					case 16: //"-tr", "--track"
						trackDependencies = parseIntegerArray(options);
						break;
					case 17: //"-tm", "--trackMask", "{string string}", "", "0",
						String[] r = parseStringArray(options, 2);
						trackMaskInput = r[0];
						trackMaskOutput = r[1];
						if(trackMaskInput.endsWith(".raw")){
							argsSpecification[18/*-mg*/][4] = "1";
						}
						break;
					case 18: //"-mg", "--maskGeometry", "{int int int int boolean}", "", "0",
						trackMaskGeometry = parseIntegerArray(options, 5);
						break;
					case 19: //"-t", "--transform", "", "0", "0",
						transform = parseIntegerPositive(options);
						break;
					}
				}else{
					throw new Exception("Argument \"" + args[argNum] + "\" repeated.");
				}
			}else{
				throw new Exception("Argument \"" + args[argNum] + "\" unrecognized.");
			}
		}

		//Check mandatory arguments
		for(int i = 0; i < argsSpecification.length; i++){
			if(argsSpecification[i][4].compareTo("1") == 0){
				if(!argsFound[i]){
					throw new Exception("Argument \"" + argsSpecification[i][0] + "\" is mandatory (\"-h\" displays help).");
				}
			}
		}
	}

	/**
	 * Finds the argument string in arguments specification array.
	 *
	 * @param arg argument to find out in argsSpecification
	 * @return the argument index of argsSpecification (-1 if it doesn't exist)
	 */
	int argFind(String arg){
		int argFound = 0;
		boolean found = false;

		while((argFound < argsSpecification.length) && !found){
			if((arg.compareTo(argsSpecification[argFound][0]) == 0) || (arg.compareTo(argsSpecification[argFound][1]) == 0)){
				found = true;
			}else{
				argFound++;
			}
		}
		return(found ? argFound: -1);
	}

	/**
	 * This function shows arguments information to console.
	 */
	public void showArgsInfo(){
		System.out.println("Arguments specification: ");
		for(int numArg = 0; numArg < argsSpecification.length; numArg++){
			char beginMandatory = '{', endMandatory = '}';
			if(argsSpecification[numArg][4].compareTo("0") == 0){
				//No mandatory argument
				beginMandatory = '[';
				endMandatory = ']';
			}
			System.out.print("\n" + beginMandatory + " ");
			System.out.print("{" + argsSpecification[numArg][0] + "|" + argsSpecification[numArg][1] + "} " + argsSpecification[numArg][2]);
			System.out.println(" " + endMandatory);
			System.out.println("  Explanation:\n    " + argsSpecification[numArg][5]);
			System.out.println("  Default value: " + argsSpecification[numArg][3]);
		}
	}


	/////////////////////
	//PARSING FUNCTIONS//
	/////////////////////
	//These functions receives a string array that contains in first position the argument and then their options//

	int parseIntegerPositive(String[] options) throws Exception{
		int value = 0;

		if(options.length == 2){
			try{
				value = Integer.parseInt(options[1]);
				if(value < 0){
					throw new Exception("\"" + options[1] + "\" of argument \"" + options[0] + "\" is must be a positive integer.");
				}
			}catch(NumberFormatException e){
				throw new Exception("\"" + options[1] + "\" of argument \"" + options[0] + "\" is not a parsable integer.");
			}
		}else{
			throw new Exception("Argument \"" + options[0] + "\" takes one option. Try \"-h\" to display help.");
		}
		return(value);
	}

	String parseString(String[] options) throws Exception{
		String value = "";

		if(options.length == 2){
			value = options[1];
		}else{
			throw new Exception("Argument \"" + options[0] + "\" takes one option. Try \"-h\" to display help.");
		}
		return(value);
	}

	String[] parseStringArray(String[] options, int numOptions) throws Exception{
		String[] value = new String[numOptions];

		if(options.length == numOptions + 1){
			for (int i = 0; i < numOptions; i++) {
				value[i] = options[1 + i];
			}
		}else{
			throw new Exception("Argument \"" + options[0] + "\" takes " + numOptions + " option. Try \"-h\" to display help.");
		}
		return(value);
	}
	
	int[] parseIntegerArray(String[] options) throws Exception{
		int[] value = null;

		if(options.length >= 2){
			value = new int[options.length - 1];
			for(int numOption = 1; numOption < options.length; numOption++){
				try{
						value[numOption - 1] = Integer.parseInt(options[numOption]);
				}catch(NumberFormatException e){
					throw new Exception("\"" + options[numOption] + "\" of argument \"" + options[0] + "\" is not a parsable integer.");
				}
			}
		}else{
			throw new Exception("Argument \"" + options[0] + "\" takes one or more options. Try \"-h\" to display help.");
		}
		return(value);
	}

	int[] parseIntegerArray(String[] options, int numOptions) throws Exception{
		int[] value = null;

		if(options.length == numOptions+1){
			value = new int[options.length - 1];
			for(int numOption = 1; numOption < options.length; numOption++){
				try{
						value[numOption - 1] = Integer.parseInt(options[numOption]);
				}catch(NumberFormatException e){
					throw new Exception("\"" + options[numOption] + "\" of argument \"" + options[0] + "\" is not a parsable integer.");
				}
			}
		}else{
			throw new Exception("Argument \"" + options[0] + "\" takes " + numOptions +" options. Try \"-h\" to display help.");
		}
		return(value);
	}
	
	float[] parseFloatArray(String[] options) throws Exception{
		return parseFloatArray(options, options.length - 1);
	}
	
	float[] parseFloatArray(String[] options, int numOptions) throws Exception{
		float[] value = null;

		if(options.length == numOptions+1){
			value = new float[options.length - 1];
			for(int numOption = 1; numOption < options.length; numOption++){
				try{
						value[numOption - 1] = Float.parseFloat(options[numOption]);
				}catch(NumberFormatException e){
					throw new Exception("\"" + options[numOption] + "\" of argument \"" + options[0] + "\" is not a parsable float.");
				}
			}
		}else{
			throw new Exception("Argument \"" + options[0] + "\" takes one or more options. Try \"-h\" to display help.");
		}
		return(value);
	}


	///////////////////////////
	//ARGUMENTS GET FUNCTIONS//
	///////////////////////////

	public String getCommandLine() {
		return(commandLine);
	}
	public String getInputImageFile(){
		return(inputImageFile);
	}
	public String getTransformFile(){
		return(transformFile);
	}
	public String getOutputImageFile(){
		return(outputImageFile);
	}
	public int[] getInputImageGeometry(){
		return(inputImageGeometry);
	}
	public int[] getOutputImageGeometry(){
		return(outputImageGeometry);
	}
	public int[] getCoefficientsApproximation(){
		return(coefficientsApproximation);
	}
	public float[] getMinValue(){
		return(minValue);
	}
	public float[] getMaxValue(){
		return(maxValue);
	}
	public int getDimension(){
		return(dimension);
	}
	public int getDirection(){
		return(direction);
	}
	public int getzBegin(){
		return(zBegin);
	}
	public boolean getEnableClustering(){
		return(enableClustering);
	}
	public int[] getClusteringOptions(){
		return(clusteringOptions);
	}
	public String getDumpEquivalentMatrix(){
		return(dumpEquivalentMatrix);
	}
	public float getEnableSubsampling(){
		return(subsampling);
	}
	public int getTransformFilter() {
		return(transform);
	}
	public int[] getTrackDependencies() {
		return trackDependencies;
	}
	public int[] getTrackMaskGeometry() {
		return trackMaskGeometry;
	}
	public String getTrackMaskInput() {
		return trackMaskInput;
	}
	public String getTrackMaskOutput() {
		return trackMaskOutput;
	}
}

