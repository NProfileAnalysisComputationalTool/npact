import java.io.*;
import java.net.URLDecoder;
import java.util.*;

public class NpactCalc {	
	private String npactDirectory;
	private File gbkFile;
	private String fileName;
	private String date;
	private String newDirectory;
	private Date resultfolder;
	
	//Extract.c keys
	private String extractCommand = new String("analysis/extract");
	private String geneDescriptorKey1 = new String("gene");
	private String geneDescriptorKey2 = new String("locus_tag");
	private String geneDescriptorSkip1 = new String("0");
	private String geneDescriptorSkip2 = new String("0");
	
	//Nprofiles.c keys
	private String nprofileCommand = new String("analysis/nprofile");
	private String window_size = new String("201");
	private String step = new String("51");
	private String period = new String("3");
	private String nucleotides = new String("cg");
	private String length;
	
	//Acgt_gamma.c keys
	private String acgtCommand = new String("analysis/acgt_gamma");
	private String run_prediction = new String("True");
	private String significance = new String("0.001");
	
	//Allplots keys
	////////////////////////////////////Array positions///////////////////////////////////
	///////////////0 - File of unbiased CDs                        ///////////////////////
	///////////////1 - File of conserved CDs                       ///////////////////////
	///////////////2 - File of new CDs                             ///////////////////////
	///////////////3 - File of published rejected CDs              ///////////////////////
	///////////////4 - File of stretches where CD is asymmetric    ///////////////////////
	///////////////5 - File of published accepted CDs              ///////////////////////
	///////////////6 - File of potential new CDs                   ///////////////////////
	///////////////7 - File of blocks from new ORFs as cds         ///////////////////////
	///////////////8 - File of blocks from annotated genes as cds  ///////////////////////
	///////////////9 - File of GeneMark regions                    ///////////////////////
	///////////////10 - File of C+G coding potential regions       ///////////////////////
	///////////////11 - File of met positions                      ///////////////////////
	///////////////12 - File of stop positions                     ///////////////////////
	///////////////13 - File of tatabox positions                  ///////////////////////
	///////////////14 - File of capbox positions                   ///////////////////////
	///////////////15 - File of ccaatbox positions                 ///////////////////////
	///////////////16 - File of gcbox positions                    ///////////////////////
	///////////////17 - File of kozak positions                    ///////////////////////
	///////////////18 - File of palindrome positions and size      ///////////////////////
	///////////////19 - File list of nucleotides in 200bp windows  ///////////////////////
	///////////////20 - File list of nucleotides in 100bp windows  ///////////////////////
	//////////////////////////////////////////////////////////////////////////////////////
	private String[] apKeys = new String[21];
	private String description = new String("None");
	private String following_page_title = new String("Page ");
	private String bp_per_page = new String("50000");
	private int start_base;
	private int end_base;
	private int total_pages;
	private String allplotsCommand = new String("analysis/Allplots");
	
	public NpactCalc()
	{
		//Initialize the Allplots key array
		for (int i = 0; i < apKeys.length; i++)
		{
			apKeys[i] = "None";
		}
		
		//Initialize the various c code command calls
		try
		{
			npactDirectory = NpactCalc.class.getProtectionDomain().getCodeSource().getLocation().getPath();
			npactDirectory = URLDecoder.decode(npactDirectory, "UTF-8");
			
			//Need to remove starting "/" if on windows
			if (System.getProperty("os.name").toLowerCase().indexOf("win") >= 0)
			{
				npactDirectory = npactDirectory.substring(1, npactDirectory.lastIndexOf("/") + 1);
			}
			else
			{
				npactDirectory = npactDirectory.substring(0, npactDirectory.lastIndexOf("/") + 1);
		
			}
		}
		catch(Exception e)
		{}
	}

	public boolean verifyGBK(File newFile)
	{	
		gbkFile = newFile;
		String temp = gbkFile.getName();
		fileName = temp.substring(0, (temp.length() - 4));
		
		//Validate the .gbk file and extract key information
		String extension = gbkFile.getName().substring((gbkFile.getName().length() - 4));
		if (!(extension.equals(".gbk")))
		{
			//Return false: will display an error message that states file type must be .gbk
			return false;
		}
		//The user correctly chose a .gbk file
		else
		{	
			//Extract key data from the .gbk file
			String[] gbkInfo = gbkExtract(gbkFile, fileName);
			
			//The line tags were not correct
			if (gbkInfo == null)
			{
				//Return false: will display an error message that states file type must be .gbk
				return false;
			}
			else
			{
				length = gbkInfo[0];
				date = gbkInfo[1];
				description = gbkInfo[2];
			
				//Return true: will enable the configuration options and set starting values
				return true;
			}
		}
	}
	
	//Function that extracts data from the first line of the .gbk file
	private String[] gbkExtract(File gbkFile, String fileName)
	{
		System.out.println("Beginning analysis of " + gbkFile.getName());
		
		try
		{
			//Create file reading variables
			FileInputStream fs = new FileInputStream(gbkFile);
			DataInputStream ds = new DataInputStream(fs);
			BufferedReader br = new BufferedReader(new InputStreamReader (ds));
				
			//Read the first line of the .gbk file
			String line = br.readLine();
			
			//Split the first line into several pieces of array data
			String[] line1 = line.split("\\s+");
				
			//Read the second line of the .gbk file
			line = br.readLine();
				
			//Split the second line into several pieces of array data
			String[] line2 = line.split("\\s+");
				
			//Make sure that the first line of the file has the appropriate tag
			if (!line1[0].equals("LOCUS"))
			{
				return null;
			}
				
			//Extract the definition line from .gbk file if it has one
			if (line2[0].equals("DEFINITION") && line2.length > 1)
			{
				description = line2[1];
				for (int i = 2; i < line2.length; i++)
				{
					description = description + " " + line2[i];
				}
					
				//Line1[2] = sequence length
				//Line1[7] = date
				String[] results = {line1[2], line1[7], description};
					
				return results;
			}
			else
			{
				//Line1[2] = sequence length
				//Line1[7] = date
				//Line2[1] = name of the file
				String[] results = {line1[2], line1[7], fileName};
					
				return results;
			}
		}
		catch (Exception e)
		{
			System.out.println("Error: " + e.getMessage());
			System.exit(1);
		}
		
		return null;
	}
	
	//Function that executes the extract.c code
	public int runExtract()
	{
		//Find better place for this later//
		//Finalize the total_page size
		total_pages = (((end_base - start_base) * 1) / Integer.parseInt(bp_per_page)) + 1;
		
		//Assemble the command array for the extract function call
		String[] extractCommandArray = {extractCommand, gbkFile.getAbsolutePath(), geneDescriptorSkip1,
				geneDescriptorKey1, geneDescriptorSkip2, geneDescriptorKey2};
		
		//Create a new process builder to execute the call to extract.c
		ProcessBuilder pb = new ProcessBuilder(extractCommandArray);
		pb.directory(new File(npactDirectory));
		int shellExitStatus = 1;
				
		try
		{
			//Open a shell and run the command
			Process shell = pb.start();
					
			//Create variables to read the output stream of extract.c
			InputStream shellIn = shell.getInputStream();
			InputStreamReader shellReader = new InputStreamReader(shellIn);
			BufferedReader bShellReader = new BufferedReader(shellReader);
			
			//Create variables to consume error stream extract.c
			InputStream errIn = shell.getErrorStream();
			InputStreamReader errReader = new InputStreamReader(errIn);
			BufferedReader bErrReader = new BufferedReader(errReader);
					
			//Read the output from extract.c and create an output file for the content
			writeFile(bShellReader, fileName, ".genes");
			
			//Consume the error lines inside the buffer so the process returns properly (necessary on Windows)
			consumeError(bErrReader);
					
			//Record the new file's location in the AP keys array
			apKeys[5] = npactDirectory + fileName + ".genes";
		
			//Wait for the process to terminate and print the exit status
			shellExitStatus = shell.waitFor();
			
			if (shellExitStatus == 0)
			{
				System.out.println("Extract operation complete.");
			}
			else
			{
				System.out.println("Error during extract operation.");
			}
			
			shellIn.close();
		}
				
		//Catch exceptions
		catch (IOException e)
		{
			System.out.println("Error occured while executing extract command. Error Description: " + 
					e.getMessage());
		}
				
		catch (InterruptedException e)
		{
			System.out.println("Error occured while executing extract command. Error Description: " +
					e.getMessage());
		}
		
		return shellExitStatus;
	}
	
	//Function that executes the nprofile.c code
	public int runNprofile ()
	{
		//Assemble the command array for the nprofile function call
		String[] nprofileCommandArray = {nprofileCommand, "-b", nucleotides, gbkFile.getAbsolutePath(),
				"1", length, window_size, step, period};
		
		//Create a new process builder to execute the call to nprofile.c
		ProcessBuilder pb = new ProcessBuilder(nprofileCommandArray);
		pb.directory(new File(npactDirectory));
		int shellExitStatus = 1;
				
		try
		{
			//Open a shell and run the command
			Process shell = pb.start();
			
			//Create variables to read the output stream of nprofile.c
			InputStream shellIn = shell.getInputStream();
			InputStreamReader shellReader = new InputStreamReader(shellIn);
			BufferedReader bShellReader = new BufferedReader(shellReader);
			
			//Create variables to consume error stream nprofile.c
			InputStream errIn = shell.getErrorStream();
			InputStreamReader errReader = new InputStreamReader(errIn);
			BufferedReader bErrReader = new BufferedReader(errReader);
					
			//Read the output from extract.c and create an output file for the content
			writeFile(bShellReader, fileName, ".nprofile");
			
			//Consume the error lines inside the buffer so the process returns properly (necessary on Windows)
			consumeError(bErrReader);
			
			//Record the new file's location in the AP keys array
			apKeys[19] = npactDirectory + fileName + ".nprofile";
					
			//Wait for the process to terminate and print the exit status
			shellExitStatus = shell.waitFor();
			
			if (shellExitStatus == 0)
			{
				System.out.println("Nprofile operation complete.");
			}
			else
			{
				System.out.println("Error during nprofile operation.");
			}
			
			shellIn.close();
		}
				
		//Catch exceptions
		catch (IOException e)
		{
			System.out.println("Error occured while executing extract command. Error Description: " + 
					e.getMessage());
		}
				
		catch (InterruptedException e)
		{
			System.out.println("Error occured while executing extract command. Error Description: " +
					e.getMessage());
		}
		
		return shellExitStatus;
	}
	
	//Function that executes the acgt_gamma.c code
	public int runAcgt ()
	{
		//Execute acgt_gamma only if the user did not check the 'skip prediction' checkbox
		if (run_prediction.equals("True"))
		{
			//Assemble the command array for the acgt_gamma function call
			String[] acgtCommandArray = {acgtCommand, "-q", "-q",  
					gbkFile.getAbsolutePath(), significance};
			
			//Create a new process builder to execute the call to acgt_gamma.c
			ProcessBuilder pb = new ProcessBuilder(acgtCommandArray);
			pb.directory(new File(npactDirectory));
			int shellExitStatus = 1;
						
			try
			{
				//Open a shell and run the command
				Process shell = pb.start();
			
				//Create variables to read the output stream of nprofile.c
				InputStream shellIn = shell.getInputStream();
					
				//Wait for the process to terminate and print the exit status
				shellExitStatus = shell.waitFor();
				
				if (shellExitStatus == 0)
				{
					System.out.println("Acgt_gamma operation complete.");
				}
				else
				{
					System.out.println("Error during acgt_gamma operation.");
				}
				
				//Record the new file locations in the AP keys array
				apKeys[2] = npactDirectory + fileName + ".newcds";
				apKeys[3] = npactDirectory + fileName + ".modified";
				apKeys[10] = npactDirectory + fileName + ".profiles";
				
				shellIn.close();
			}
						
			//Catch exceptions
			catch (IOException e)
			{
				System.out.println("Error occured while executing extract command. Error Description: " + 
						e.getMessage());
			}
						
			catch (InterruptedException e)
			{
				System.out.println("Error occured while executing extract command. Error Description: " +
						e.getMessage());
			}
		
			return shellExitStatus;
		}
		
		return 0;
	}
	
	//Function that executes Allplots.c
	public int[] runAllplots(int page_num)
	{
		//Declarations
		int results[] = {2, total_pages};
		
		//If this holds true, we need to keep running allplots
		if (start_base < end_base)
		{
			//Create the .def file for the current .ps page
			writeAllplots(page_num);
			
			System.out.println("Generating page " + page_num + " out of " + total_pages);
			
			//Create a new process builder to execute the call to Allplots.c
			String[] command1 = {allplotsCommand, Integer.toString(start_base), bp_per_page, "5", "1000", period, Integer.toString(end_base)};
			String[] command2 = {allplotsCommand, Integer.toString(start_base), bp_per_page, "5", "1000", period};
			ProcessBuilder pb;
			
			if (start_base + Integer.parseInt(bp_per_page) >= end_base)
			{
				pb = new ProcessBuilder(command1);
			}
			else
			{
				pb = new ProcessBuilder(command2);
			}
			pb.directory(new File(npactDirectory));
			
			try
			{
				//Open a shell and run the command
				Process shell = pb.start();
				
				//Create variables to read the output stream of Allplots.c
				InputStream shellIn = shell.getInputStream();
				InputStreamReader shellReader = new InputStreamReader(shellIn);
				BufferedReader bShellReader = new BufferedReader(shellReader);
				
				//Create variables to consume error stream Allplots.c
				InputStream errIn = shell.getErrorStream();
				InputStreamReader errReader = new InputStreamReader(errIn);
				BufferedReader bErrReader = new BufferedReader(errReader);
						
				//Read the output from extract.c and create an output file for the content
				writeFile(bShellReader, ".result" + page_num, ".ps");
				
				//Consume the error lines inside the buffer so the process returns properly (necessary on Windows)
				consumeError(bErrReader);
						
				//Wait for the process to terminate and print the exit status
				int shellExitStatus = shell.waitFor();				
				
				if (shellExitStatus != 0)
				{
					System.out.println("Error during allplots operation.");
				}
				
				
				shellIn.close();
				shellReader.close();
				bShellReader.close();
				results[0] = shellExitStatus;
			}
									
			//Catch exceptions
			catch (IOException e)
			{
				System.out.println("Error occured while executing allplots command. Error Description: " + 
						e.getMessage());
			}
									
			catch (InterruptedException e)
			{
				System.out.println("Error occured while executing allplots command. Error Description: " +
						e.getMessage());
			}
			
			start_base += Integer.parseInt(bp_per_page);
		}
		
		return results;
	}
	
	//Function that writes the Allplots.def file
	private void writeAllplots(int page_num)
	{
		//Create the Allplots.def file
		try
		{
			FileWriter fstream = new FileWriter(npactDirectory + "Allplots.def");
			BufferedWriter out = new BufferedWriter(fstream);
	
			//Write the first four lines of the file
			out.write("FOOBAR " + length + "\n"); //Plot title
			for (int i = 0; i < nucleotides.length(); i++) //Nucleotide(s) plotted
			{
				if (i != 0)
				{
					out.write(" + ");
				}
				out.write(nucleotides.charAt(i));
			}
			out.write("\n");
			out.write(description + "\n");
			out.write(following_page_title + page_num + "\n");
			
			//Print out the AP key array to the file
			for (int i = 0; i < apKeys.length; i++)
			{
				out.write(apKeys[i] + "\n");
			}
			
			out.close();
		}
		
		catch (IOException e)
		{
			System.out.println("An error occured while writing the Allplots.def file: " +
					e.getMessage());
		}
	}
	
	//Function to combine the various PS output files from Allplots.c
	public void combineOutput()
	{	
		System.out.println("Combining output files into single file.");
		
		//Create the new output file
		try 
		{
			FileWriter fstream = new FileWriter(npactDirectory + 
					gbkFile.getName().substring(0, (gbkFile.getName().length() - 4)) + ".ps");
			BufferedWriter out = new BufferedWriter(fstream);
			
			//Write the header lines for the output file
			out.write("%!PS-Adobe-2.0\n\n");
			out.write("%%Pages: " + total_pages +"\n\n");
			out.write("%%EndComments" + "\n");
			
			//Loop through all of the numbered output files, placing their output into the new
			//combined output
			for (int i = 1; i <= total_pages; i++)
			{
				//Create file reading variables
				FileInputStream fs = new FileInputStream(npactDirectory + ".result" + i + ".ps");
				DataInputStream ds = new DataInputStream(fs);
				BufferedReader br = new BufferedReader(new InputStreamReader (ds));
				String line;
				
				//Clear out the first two lines of the selected file
				br.readLine();
				br.readLine();
				
				//Write out line indicating what file we are on
				out.write("%%Page: " + i + "\n");
				
				while ((line = br.readLine()) != null)
				{
					out.write(line + "\n");
				}
				
				if (i == total_pages)
				{
					//Write end of file line
					out.write("%%Trailer" + "\n");
					out.write("%%Pages: " + total_pages + "\n");
					out.write("%%EOF");
				}
				
				//Close the streams
				fs.close();
				ds.close();
				br.close();
			}
			
			//Close the writing stream
			out.close();
		}
		
		catch (IOException e)
		{
			System.out.println("An error occured while writing the combined results file: " +
					e.getMessage());
		}
	}
	
	//Function used to delete all extraneous files and relocate results files into proper final directory
	public void cleanOutput()
	{		
		//Loop through all of the fileNum number of intermediate output files and delete them
		for (int i =1; i <= total_pages; i++)
		{
			//Access the current result file and delete it
			File temp = new File(npactDirectory + ".result" + i + ".ps");
			temp.delete();
		}
		
		//Delete the Allplots.def file
		File temp = new File(npactDirectory + "Allplots.def");
		temp.delete();
		
		//Get the current date/time to use as the name of a new directory
		resultfolder = new Date();
		
		//If the user ran the program from the command line, newDirectory is null and needs to be set to the current directory
		if (newDirectory == null)
		{
			newDirectory = System.getProperty("user.dir");
		}
		
		//Change the direction of the slash depending on what operating system we are on
		//Replace the colons in the resultfolder string so Windows can create the folders
		if (System.getProperty("os.name").toLowerCase().indexOf("win") >= 0)
		{
			newDirectory = newDirectory + "\\" + resultfolder.toString().replace(':', '.');
		}
		else
		{
			newDirectory = newDirectory + "/" + resultfolder.toString().replace(':', '.');
		}
		
		//Create a new directory and move all of the result files into it
		try
		{
			if (new File(newDirectory).mkdir())
			{
				System.out.println("Results directory created: " + newDirectory + "\n");
				String oldDirectory = new String(npactDirectory);
			
				//Add the appropriate type of slash at the end of newDirectory depending on what type of OS we are on
				if (System.getProperty("os.name").toLowerCase().indexOf("win") >= 0)
				{
					newDirectory = newDirectory + "\\";
				}
				else
				{
					newDirectory = newDirectory + "/";
				}
				
				File altcds = new File(oldDirectory + fileName + ".altcds");
				altcds.renameTo(new File(newDirectory + altcds.getName()));
				
				File confirmed = new File(oldDirectory + fileName + ".confirmed");
				confirmed.renameTo(new File(newDirectory + confirmed.getName()));
				
				File excluded = new File(oldDirectory + fileName + ".excluded");
				excluded.renameTo(new File(newDirectory + excluded.getName()));
				
				File modified = new File(oldDirectory + fileName + ".modified");
				modified.renameTo(new File(newDirectory + modified.getName()));
				
				File newcds = new File(oldDirectory + fileName + ".newcds");
				newcds.renameTo(new File(newDirectory + newcds.getName()));
				
				File profiles = new File(oldDirectory + fileName + ".profiles");
				profiles.renameTo(new File(newDirectory + profiles.getName()));
				
				File repetitive = new File(oldDirectory + fileName + ".repetitive");
				repetitive.renameTo(new File(newDirectory + repetitive.getName()));
				
				File genes = new File(oldDirectory + fileName + ".genes");
				genes.renameTo(new File(newDirectory + genes.getName()));
				
				File nprofile = new File(oldDirectory + fileName + ".nprofile");
				nprofile.renameTo(new File(newDirectory + nprofile.getName()));
				
				File verbose = new File(oldDirectory + fileName + ".verbose");
				verbose.renameTo(new File(newDirectory + verbose.getName()));
				
				File results = new File(oldDirectory + gbkFile.getName().substring(0, (gbkFile.getName().length() - 4)) +
						".ps");
				results.renameTo(new File(newDirectory + results.getName()));
				
				newDirectory = newDirectory + results.getName();
			}
		}
		
		catch(Exception e)
		{
			System.out.println("Error attempting to create output file directory: " + e.getMessage());
		}
	}
	
	//Function that creates output files from the output of the various c files
	private void writeFile(BufferedReader outReader, String fileName, String extension)
	{
		String line;
		
		try
		{
			//Create the output file
			FileWriter fstream = new FileWriter(npactDirectory + fileName + extension);
			BufferedWriter out = new BufferedWriter(fstream);
				
			//Write the output to the file
			while((line = outReader.readLine()) != null)
			{
				out.write(line + "\n");
			}
				
			//Close the output stream
			out.close();
		}
			
		catch (IOException e)
		{
			System.out.println("Error occured while executing extract command. Error Description: " + 
					e.getMessage());
		}
	}
	
	//Function that reads and consumes error lines printed as output from the various c files
	private void consumeError(BufferedReader errReader)
	{	
		try
		{
			while (errReader.readLine() != null)
			{
			
			}
		}
		catch (Exception e)
		{
			System.out.println("Error while consuming c code error output.");
		}
	}
	
	//Function used to delete all output files should an error have occurred during execution
	public void deleteOutput(int allplotsFlag)
	{
		String directory = new String(npactDirectory);
		
		File temp = new File(directory + fileName + ".altcds");
		temp.delete();
			
		temp = new File(directory + fileName + ".confirmed");
		temp.delete();
			
		temp = new File(directory + fileName + ".excluded");
		temp.delete();
			
		temp = new File(directory + fileName + ".modified");
		temp.delete();
			
		temp = new File(directory + fileName + ".newcds");
		temp.delete();
			
		temp = new File(directory + fileName + ".profiles");
		temp.delete();
			
		temp = new File(directory + fileName + ".repetitive");
		temp.delete();
			
		temp = new File(directory + fileName + ".genes");
		temp.delete();
			
		temp = new File(directory + fileName + ".nprofile");
		temp.delete();
			
		temp = new File(directory + fileName + ".verbose");
		temp.delete();
			
		//If allplots.c made results files, we need to delete all of those them
		if (allplotsFlag > 0)
		{
			for (int i = 1; i <= allplotsFlag; i++)
			{
				temp = new File(directory + "result" + i + ".ps");
				temp.delete();
			}
		}
	}
	
	//////////////////////
	//Get & Set methods //
	//////////////////////
	public void resetGBKFile()
	{
		this.gbkFile = null;
	}
	
	public String getNewDirectory()
	{
		return this.newDirectory;
	}
	
	public void setNewDirectory(String text)
	{
		this.newDirectory = text;
	}
	
	public Date getResultFolder()
	{
		return this.resultfolder;
	}
	
	public String getDescription()
	{
		return description;
	}
	
	public void setDescription(String text)
	{
		this.description = text;
	}
	
	public String getFollowingTitle()
	{
		return following_page_title;
	}
	
	public void setFollowingTitle(String text)
	{
		this.following_page_title = text;
	}
	
	public String getLength()
	{
		return length;
	}
	
	public void setLength(String text)
	{
		this.length = text;
	}
	
	public String getNucleotides()
	{
		return nucleotides;
	}
	
	public void nucleotideReset()
	{
		this.nucleotides = "";
	}
	
	public void nucleotideConcat(String text)
	{
		this.nucleotides = this.nucleotides + text;
	}
	
	public String getSkip()
	{
		return run_prediction;
	}
	
	public void setSkip(String text)
	{
		this.run_prediction = text;
	}
	
	public String getSignificance()
	{
		return significance;
	}
	
	public void setSignificance(String text)
	{
		this.significance = text;
	}
	
	public int getStartBase()
	{
		return start_base;
	}
	
	public void setStartBase(int value)
	{
		this.start_base = value;
	}
	
	public int getEndBase()
	{
		return end_base;
	}
	
	public void setEndBase(int value)
	{
		this.end_base = value;
	}
}
