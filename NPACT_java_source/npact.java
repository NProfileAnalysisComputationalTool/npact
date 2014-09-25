import java.io.*;

public class npact {
	//Initialize the calculation object that will handle all processing requests
	static NpactCalc calculation = new NpactCalc();
	
	public static void main(String args[])
	{	
		//If the user didn't pass in any arguments on the command line, start the GUI version of the program
		if (args.length == 0)
		{
			NpactGUI gui = new NpactGUI();
			gui.launchFrame();
		}
		//Run command line version of NPACT if command line arguments were passed in
		else
		{
			run_commandline(args);
		}
	}
	
	//Function that handles command-line execution of the program
	public static void run_commandline(String[] args)
	{		
		//If the first argument is the help delimiter, print out the help screen and exit
		if (args[0].equals("-h"))
		{
			//Print out the help screen and exit
			printHelp();
			System.exit(0);
		}
		
		//Verify that the first argument passed in is a valid file
		File gbkFile = new File(args[0]);
		
		if (!gbkFile.isFile())
		{
			System.out.println("A file with that name does not exist.");
			System.exit(1);
		}
		
		//Verify that the first argument passed in is a .gbk file
		if (calculation.verifyGBK(gbkFile))
		{
			//The user supplied a .gbk file, continue with process
			//Process all of the  options the user inputed
			validateInput(args);
			
			//Used to keep track of the exit status of the various c code calls
			int exitStatus = 0;
			
			//Run the extract.c file and delete results if there was an error
			exitStatus = calculation.runExtract();
			if (exitStatus != 0)
			{
				//Delete all files that may have been made during execution
				calculation.deleteOutput(0);
				System.out.println("There was an error during the extract phase of execution.");
				System.exit(1);
			}
			
			//Call the function that executes the nprofile command
			exitStatus = calculation.runNprofile();
			if (exitStatus != 0)
			{
				//Delete all files that may have been made during execution
				calculation.deleteOutput(0);
				System.out.println("There was an error during the nprofile phase of execution.");
				System.exit(1);
			}
			
			//Call the function that executes the acgt_gamma command
			calculation.runAcgt();
			if (exitStatus != 0)
			{
				//Delete all files that may have been made during execution
				calculation.deleteOutput(0);
				System.out.println("There was an error during the prediction phase of execution.");
				System.exit(1);
			}
			
			//Execute Allplots.c
			int[] allplotsResults = {0, 0};
			int page_num = 1;
			
			while (allplotsResults[0] < 2)
			{
				allplotsResults = calculation.runAllplots(page_num);
				page_num++;
				
				if (allplotsResults[0] == 1)
				{
					//Delete all files that may have been made during execution
					calculation.deleteOutput(allplotsResults[1]);
					
					//Display an error message
					System.out.println("There was an error during the plotting phase of execution.");
					System.exit(1);
				}
			}
			
			//Combine the results PS files into a single file
			calculation.combineOutput();
			
			//Delete the files used to make the combined PS file
			//Relocate the used files to an output directory
			calculation.cleanOutput();
		}
		else
		{
			//Display an error message that states file type must be .gbk
			System.out.println("The file you have chosen is not a .gbk file.");
			System.exit(1);
		}
	}
	
	//Function that validates the input string the user supplied to the program
	public static void validateInput(String[] args)
	{
		//Set the start_base and end_base values as they are not set naturally at start
		calculation.setStartBase(0);
		calculation.setEndBase(Integer.parseInt(calculation.getLength()));
		
		String flag = new String();
		
		//Iterate over all of the args passed into the program, skipping the .gbk file
		for (int i = 1; i < args.length; i++)
		{
			//If the i is odd then an option flag needs to be in this position
			if (i % 2 == 1)
			{
				flag = args[i];
				
				if (flag.equals("-h"))
				{
					printHelp();
					System.exit(0);
				}
			}
			//If the i is even then we need to set whatever option flag has been passed to this value
			else
			{
				//First-page flag is set
				if (flag.equals("-firstp"))
				{
					calculation.setDescription(args[i]);
				}
				//Following-page flag is set
				else if (flag.equals("-followp"))
				{
					calculation.setFollowingTitle(args[i]);
				}
				//Length flag is set
				else if (flag.equals("-l"))
				{
					calculation.setLength(args[i]);
				}
				//Nucleotides flag is set
				else if (flag.equals("-n"))
				{
					//Verify that the nucleotides were entered in a correct fashion
					validateNucleotides(args[i]);
				}
				//Skip prediction flag is set
				else if (flag.equals("-sp"))
				{
					//If the user entered 't', the user wants to skip the prediction
					if (args[i].toLowerCase().equals("t"))
					{
						calculation.setSkip("False");
					}
					//If the user entered 'f', the user wants to run the prediction
					else if (args[i].toLowerCase().equals("f"))
					{
						calculation.setSkip("True");
					}
					//If the user entered anything else, quit
					else
					{
						System.out.println("Please only use 't' for true or 'f' for false when setting the skip prediciton flag.");
						System.exit(1);
					}
				}
				//Prediction significance flag is set
				else if (flag.equals("-ps"))
				{
					//Continue execution if the user entered one of the three available prediction significances,
					//otherwise exit
					if (Float.parseFloat(args[i]) == 0.01)
					{
						calculation.setSignificance("0.01");
					}
					else if (Float.parseFloat(args[i]) == 0.001)
					{
						calculation.setSignificance("0.001");
					}
					else if (Float.parseFloat(args[i]) == 0.0001)
					{
						calculation.setSignificance("0.0001");
					}
					else
					{
						System.out.println("Invalid significance value entered. Please enter either 0.01, 0.001, or 0.0001");
						System.exit(1);
					}
				}
				//Start base flag is set
				else if (flag.equals("-sb"))
				{
					calculation.setStartBase(Integer.parseInt(args[i]));
				}
				//End base flag is set
				else if (flag.equals("-eb"))
				{
					calculation.setEndBase(Integer.parseInt(args[i]));
				}
				//If none of the above values for the flag are set, print help, and exit as the user has entered an 
				//invalid flag value
				else
				{
					printHelp();
					System.exit(1);
				}
			}
		}
	}
	
	//Function that validates whether the user supplied a valid nucleotide string
	public static void validateNucleotides(String nucleotides)
	{
		//Loop through the string making sure all entries are a,c,g, or t
		nucleotides = nucleotides.toLowerCase();
		for (int i = 0; i < nucleotides.length(); i++)
		{
			//The user entered a character that is invalid
			if (nucleotides.charAt(i) != 'a' && nucleotides.charAt(i) != 'c' && nucleotides.charAt(i) != 'g' && 
					nucleotides.charAt(i) != 't')
			{
				//Print a message to the user and exit the program
				System.out.println("For the nucleotide option, only enter which nucleotides you would like processed. " + 
				"For example: \"cg\" or \"acgt\"");
				System.exit(1);
			}
			//The user entered a valid nucleotide, set the appropriate flag in the calculation object
			else
			{
				calculation.nucleotideConcat(Character.toString(nucleotides.charAt(i)));
			}
		}
	}
	
	public static void printHelp()
	{
		System.out.println("NPACT Command-line Usage:\n");
		System.out.println("\tjava -jar NPACT.jar {filename}\n");
		System.out.println("Optional arguments:\n");
		System.out.println("\t-firstp\t\tthe title used for the first page of the results file");
		System.out.println("\t\t\tDefault: Name of the sequence being analyzed\n");
		System.out.println("\t-followp\tthe title used for all other pages after the first");
		System.out.println("\t\t\tDefault: Page\n");
		System.out.println("\t-l\t\tthe total length of the sequence being submitted for analysis\n");
		System.out.println("\t-n\t\tthe nucleotides to be examined during analysis");
		System.out.println("\t\t\tDefault: cg\n");
		System.out.println("\t-sp\t\tshould the prediction portion of the analysis be skipped?");
		System.out.println("\t\t\tDefault: f\n");
		System.out.println("\t-ps\t\tset the prediction signifance value to either 0.01, 0.001, or 0.0001");
		System.out.println("\t\t\tDefault: 0.001\n");
		System.out.println("\t-sb\t\tthe base from which analysis starts");
		System.out.println("\t\t\tDefault: 0\n");
		System.out.println("\t-b\t\tthe base at which analysis ends");
		System.out.println("\t\t\tDefault: last base of the sequence\n");
		System.out.println("\t-h\t\tprint out the help message for command-line usage\n");
	}
}
