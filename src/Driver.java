import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;


public class Driver {
	public static void main(String[] args) throws IOException {
		markovNetworkTest(args);
		//naiveBayesTest(args);
	}
	
	public static void markovNetworkTest(String[] args) {
		Scanner scan;
		int numVariables = 0;
		int numCliques = 0;
		int cardinalities[] = {};
		int variables[] = {};
		ArrayList<Factor> factors = new ArrayList<>();
		
		// Check for proper number of arguments (1).
		if(args.length != 1) {
			System.out.println("Usage: java -jar mn-partition.jar [path/to/network/file]");
			System.exit(0);
		}
		try {
			scan = new Scanner(new File(args[0]));
			int lineNum = 0;
			// Pass MARKOV line
			scan.nextLine();
			
			// Get number of variables
			numVariables = scan.nextInt();
			
			// Initialize arrays
			variables = new int[numVariables];
			cardinalities = new int[numVariables];
			
			// Read in cardinalities
			for(int i = 0; i < numVariables; i++) {
				cardinalities[i] = scan.nextInt();
				variables[i] = i;
			}
			
			// Get number of cliques
			numCliques = scan.nextInt();
			int cliqueSize = 0;
			
			// Initialize factors with included variables
			for(int i = 0; i < numCliques; i++) {
				cliqueSize = scan.nextInt();
				int[] vars = new int[cliqueSize];
				for(int j = 0; j < cliqueSize; j++) {
					vars[j] = scan.nextInt();
				}
				factors.add(new Factor(vars, cardinalities));
			}
			
			// Populate factors with actual values.
			for(int i = 0; i < numCliques; i++) {
				int scope = scan.nextInt();
				double[] vals = new double[scope];
				for(int j = 0; j < scope; j++) {
					vals[j] = scan.nextDouble();
				}
				factors.get(i).intializePhi(scope);
				factors.get(i).addToPhi(vals);
			}

		} catch (FileNotFoundException e) {
			// Could not find given test file so quit.
			System.out.println("Specified network file not found, quitting program.");
			System.exit(0);
		}

		// Create network
		MarkovNetwork mn = new MarkovNetwork(variables, cardinalities, factors);
		
		// Run loopy belief propogation.
		mn.loopyBP();
		
		
	}

}
