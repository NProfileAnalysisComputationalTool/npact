import javax.swing.*;
import java.io.*;
import java.net.URLDecoder;

import java.awt.*;
import java.awt.event.*;

public class NpactGUI implements ActionListener {	
	//Initialize swing frames and panels
	private JFrame frame = new JFrame("NPACT");
	private JPanel containerPanel = new JPanel();
	private JPanel boxPanel = new JPanel();
	private JPanel filePanel = new JPanel();
	private JPanel optionPanel = new JPanel();
	private JPanel nucleotidePanel = new JPanel();
	private JPanel basePanel = new JPanel();
	private JPanel startPanel = new JPanel();
	private JPanel statusPanel = new JPanel();
	private JPanel directoryPanel = new JPanel();
	private JPanel returnPanel = new JPanel();
	private NpactBKG background = new NpactBKG();
	
	//Label objects
	private JLabel fileLabel = new JLabel("Select an input file:");
	private JLabel optionLabel = new JLabel("NPACT Analysis Options:");
	private JLabel firstpgLabel = new JLabel("First page title:");
	private JLabel followpgLabel = new JLabel("Following page title:");
	private JLabel lengthLabel = new JLabel("Length:");
	private JLabel nucleotideLabel = new JLabel("Nucleotides:");
	private JLabel skipLabel = new JLabel("Skip prediction:");
	private JLabel significanceLabel = new JLabel("Prediction significance:");
	private JLabel baseLabel1 = new JLabel("Start / End base:");
	private JLabel baseLabel2 = new JLabel("/");
	private JLabel statusLabel = new JLabel();
	private JLabel directoryLabel1 = new JLabel();
	
	//Button objects
	private JButton fileButton = new JButton("Choose file");
	private JButton startButton = new JButton("Start Analysis");
	private JButton resultButton = new JButton("Show Plot");
	private JButton returnButton = new JButton("New Analysis");
	
	//Other swing objects
	private JTextField fileTextField = new JTextField();
	private JTextField firstpgTextField = new JTextField();
	private JTextField followpgTextField = new JTextField();
	private JTextField lengthTextField = new JTextField();
	private JTextField sBaseTextField = new JTextField();
	private JTextField eBaseTextField = new JTextField();
	private JCheckBox aCBox = new JCheckBox("A");
	private JCheckBox cCBox = new JCheckBox("C");
	private JCheckBox gCBox = new JCheckBox("G");
	private JCheckBox tCBox = new JCheckBox("T");
	private JCheckBox skipCBox = new JCheckBox();
	private String[] signOptions = {"0.01", "0.001", "0.0001"};
	private JComboBox significanceBox = new JComboBox(signOptions);
	private JFileChooser fc = new JFileChooser();
	
	//Initialize calculation module
	NpactCalc calculation = new NpactCalc();
	
	
	//Constructor for the GUI
	public NpactGUI()
	{		
		//Setup the file selection swing objects
		fileButton.addActionListener(this);
		fileButton.setFont(new Font("Default", Font.BOLD, 11));
		fileTextField.setPreferredSize(new Dimension (350, 20));
		fileTextField.setEditable(false);
		filePanel.add(fileButton);
		filePanel.add(fileTextField);
		filePanel.setAlignmentX(Component.LEFT_ALIGNMENT);
		filePanel.setOpaque(false);

		
		//Setup the frame
		frame.setLayout(null);
		background.setBounds(0, 0, 640, 480);
		
		frame.getContentPane().add(background);
		background.setLayout(null);
		
		background.add(containerPanel);
		containerPanel.setBounds(75, 35, 490, 445);
		containerPanel.setOpaque(false);
		containerPanel.setLayout(new FlowLayout(FlowLayout.LEFT));
		containerPanel.add(boxPanel);
		
		//Add the file selection objects to the frame
		boxPanel.setOpaque(false);
		boxPanel.setLayout(new BoxLayout(boxPanel, BoxLayout.Y_AXIS));
		boxPanel.add(fileLabel);
		boxPanel.add(filePanel);
		boxPanel.add(Box.createVerticalStrut(20));
		boxPanel.add(optionLabel);
		
		//Setup and add the panel containing the configuration options to the frame
		boxPanel.add(optionPanel);
		optionPanel.setLayout(new GridLayout(7, 7, 0, 10));
		optionPanel.setBorder(BorderFactory.createEmptyBorder(10, 50, 0, 0));
		optionPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
		optionPanel.setOpaque(false);
		optionPanel.setMaximumSize(new Dimension(400,800));
		
		//Set text field options
		firstpgTextField.setEditable(false);
		followpgTextField.setEditable(false);
		lengthTextField.setEditable(false);
		
		//Setup the nucleotide panel
		aCBox.setFont(new Font("Default", Font.BOLD, 8));
		cCBox.setFont(new Font("Default", Font.BOLD, 8));
		gCBox.setFont(new Font("Default", Font.BOLD, 8));
		tCBox.setFont(new Font("Default", Font.BOLD, 8));
		aCBox.setEnabled(false);
		cCBox.setEnabled(false);
		gCBox.setEnabled(false);
		tCBox.setEnabled(false);
		
		aCBox.setPreferredSize(new Dimension (38,19));
		cCBox.setPreferredSize(new Dimension (38,19));
		gCBox.setPreferredSize(new Dimension (38,19));
		tCBox.setPreferredSize(new Dimension (38,19));
		
		/*if (System.getProperty("os.name").indexOf("Mac") >= 0)
		{
			aCBox.setPreferredSize(new Dimension (38,19));
			cCBox.setPreferredSize(new Dimension (38,19));
			gCBox.setPreferredSize(new Dimension (38,19));
			tCBox.setPreferredSize(new Dimension (38,19));
		}
		else
		{
			aCBox.setPreferredSize(new Dimension (38,19));
			cCBox.setPreferredSize(new Dimension (38,19));
			gCBox.setPreferredSize(new Dimension (38,19));
			tCBox.setPreferredSize(new Dimension (38,19));
		}*/
		
		nucleotidePanel.add(aCBox);
		nucleotidePanel.add(cCBox);
		nucleotidePanel.add(gCBox);
		nucleotidePanel.add(tCBox);
		nucleotidePanel.setOpaque(false);
		aCBox.setOpaque(false);
		cCBox.setOpaque(false);
		gCBox.setOpaque(false);
		tCBox.setOpaque(false);
		nucleotidePanel.setPreferredSize(new Dimension(100, 5));
		nucleotidePanel.setAlignmentX(Component.LEFT_ALIGNMENT);
		nucleotidePanel.setAlignmentY(Component.TOP_ALIGNMENT);
		aCBox.setAlignmentY(Component.TOP_ALIGNMENT);
		aCBox.setAlignmentX(Component.LEFT_ALIGNMENT);
		
		//Setup the skip prediction object
		skipCBox.setOpaque(false);
		skipCBox.setEnabled(false);
		skipCBox.setPreferredSize(new Dimension (30,15));
		
		//Setup the significance options object
		significanceBox.setEnabled(false);
		significanceBox.setSelectedIndex(1);
		
		//Setup the start/end base objects
		sBaseTextField.setPreferredSize(new Dimension (70, 19));
		sBaseTextField.setEditable(false);
		eBaseTextField.setPreferredSize(new Dimension (70, 19));
		eBaseTextField.setEditable(false);
		
		basePanel.add(sBaseTextField);
		basePanel.add(baseLabel2);
		basePanel.add(eBaseTextField);
		basePanel.setOpaque(false);
		basePanel.setPreferredSize(new Dimension(100, 5));
		basePanel.setAlignmentX(Component.LEFT_ALIGNMENT);
		
		//Add various objects to the options container panel
		optionPanel.add(firstpgLabel);
		optionPanel.add(firstpgTextField);
		optionPanel.add(followpgLabel);
		optionPanel.add(followpgTextField);
		optionPanel.add(lengthLabel);
		optionPanel.add(lengthTextField);
		optionPanel.add(nucleotideLabel);
		optionPanel.add(nucleotidePanel);
		optionPanel.add(skipLabel);
		optionPanel.add(skipCBox);
		optionPanel.add(significanceLabel);
		optionPanel.add(significanceBox);
		optionPanel.add(baseLabel1);
		optionPanel.add(basePanel);
		
		//Change the label text to be white
		fileLabel.setForeground(Color.WHITE);
		optionLabel.setForeground(Color.WHITE);
		firstpgLabel.setForeground(Color.WHITE);
		followpgLabel.setForeground(Color.WHITE);
		lengthLabel.setForeground(Color.WHITE);
		nucleotideLabel.setForeground(Color.WHITE);
		skipLabel.setForeground(Color.WHITE);
		significanceLabel.setForeground(Color.WHITE);
		baseLabel1.setForeground(Color.WHITE);
		aCBox.setForeground(Color.WHITE);
		cCBox.setForeground(Color.WHITE);
		gCBox.setForeground(Color.WHITE);
		tCBox.setForeground(Color.WHITE);
		baseLabel2.setForeground(Color.WHITE);
		statusLabel.setForeground(Color.WHITE);
		directoryLabel1.setForeground(Color.WHITE);
		
		//Set the tooltip text of the options container objects
		firstpgTextField.setToolTipText("The title of the page containing the beginning of the genome.");
		followpgTextField.setToolTipText("The title of the pages after the first. Text entered here will be followed by the page number");
		lengthTextField.setToolTipText("The length, in base pairs, of the genome being analyzed.");
		aCBox.setToolTipText("The bases to count the frequency of on the primary strand.");
		cCBox.setToolTipText("The bases to count the frequency of on the primary strand.");
		gCBox.setToolTipText("The bases to count the frequency of on the primary strand.");
		tCBox.setToolTipText("The bases to count the frequency of on the primary strand.");
		skipCBox.setToolTipText("Should the acgt_gamma prediction be run? (click here to skip)");
		significanceBox.setToolTipText("What should the acgt_gamma prediction consider significant?");
		sBaseTextField.setToolTipText("The base pair coordinate at which to start graphing.");
		eBaseTextField.setToolTipText("The base pair coordinate at which to end graphing.");
		
		//Add the analysis button to the frame
		startPanel.setOpaque(false);
		startPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
		startPanel.add(startButton);
		startButton.setEnabled(false);
		startButton.addActionListener(this);
		
		boxPanel.add(Box.createVerticalStrut(30));
		boxPanel.add(startPanel);
		
		//Setup the status screen objects
		statusPanel.add(statusLabel);
		statusPanel.setOpaque(false);
		directoryPanel.add(directoryLabel1);
		directoryPanel.setOpaque(false);
		returnPanel.add(resultButton);
		returnPanel.add(returnButton);
		returnPanel.setOpaque(false);
		resultButton.addActionListener(this);
		returnButton.addActionListener(this);
		statusLabel.setIcon(new ImageIcon(new ImageIcon(NpactGUI.class.getResource("img/spinner.gif")).getImage().getScaledInstance(25, 25, Image.SCALE_DEFAULT)));
		statusLabel.setText("Analysis is being performed");
		screenSwitch("working");
		screenSwitch("menu");
	}
	
	public void launchFrame()
	{
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.setSize(640, 480);
		frame.setResizable(false);
		frame.setVisible(true);
	}
	
	public void actionPerformed(ActionEvent event)
	{
		//Handle "file" button press
		if (event.getSource() == fileButton)
		{	
			//Open a file-selection dialog box
			int returnValue = fc.showOpenDialog(this.frame);
			
			//If a file was chosen
			if (returnValue == JFileChooser.APPROVE_OPTION)
			{
				//Verify whether the file provided was an appropriate .gbk file
				if (calculation.verifyGBK(fc.getSelectedFile()))
				{
					//Display the path to the .gbk file in the text field
					fileTextField.setText(fc.getSelectedFile().getPath());
					
					//Enable the configuration options and set starting values
					enableConfig();
				}
				else
				{
					//Display an error message that states file type must be .gbk
					JOptionPane.showMessageDialog(frame, "The file you have chosen is not a .gbk file.",
							"Filetype Error", JOptionPane.ERROR_MESSAGE);
				}
			}
		}
			
		//Handle "start analysis" button press
		else if (event.getSource() == startButton)
		{
			//Change what kind of dialog is made depending on what OS is running
			if (System.getProperty("os.name").indexOf("Mac") >= 0)
			{
				//Open a Mac native file chooser to allow for directory creation
				System.setProperty("apple.awt.fileDialogForDirectories", "true");
				FileDialog d = new FileDialog(frame, "Save Results To");
			    d.setVisible(true);
			    
			    File directory = new File(d.getDirectory() + d.getFile());
			    System.setProperty("apple.awt.fileDialogForDirectories", "false");
			    
			    if (directory.exists())
			    {
			    	calculation.setNewDirectory(directory.getAbsolutePath());
			    }
			    else
			    {
			    	return;
			    }
			}
			else
			{
				//Open up a file chooser dialog box so the user can choose where to save the results
				JFileChooser saveChooser = new JFileChooser();
				saveChooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
				saveChooser.setDialogTitle("Save Results To");
				
				if(saveChooser.showOpenDialog(frame) == JFileChooser.APPROVE_OPTION )
			    {
			        calculation.setNewDirectory(saveChooser.getSelectedFile().getAbsolutePath());
			    }
				else
				{
					return;
				}
			}
			
			//Update the appropriate variables with any configuration changes made by the user
			calculation.setDescription(firstpgTextField.getText());
			calculation.setFollowingTitle(followpgTextField.getText());
			calculation.setLength(lengthTextField.getText());
			calculation.nucleotideReset();
			if (aCBox.isSelected())
			{
				calculation.nucleotideConcat("a");
			}
			if (cCBox.isSelected())
			{
				calculation.nucleotideConcat("c");
			}
			if (gCBox.isSelected())
			{
				calculation.nucleotideConcat("g");
			}
			if (tCBox.isSelected())
			{
				calculation.nucleotideConcat("t");
			}
			if (skipCBox.isSelected())
			{
				calculation.setSkip("False");
			}
			calculation.setSignificance((String)significanceBox.getSelectedItem());
			calculation.setStartBase(Integer.parseInt(sBaseTextField.getText()));
			calculation.setEndBase(Integer.parseInt(eBaseTextField.getText()));
				
			//Switch the frame over to use the status objects
			screenSwitch("working");
				
			//Create a SwingWorker instance to execute the initial c code calls without interrupting the GUI
			SwingWorker<Void, Void> task = new SwingWorker<Void, Void>() {
				public Void doInBackground() {
					//Used to keep track of the exit status of the various C code calls
					int exitStatus = 0;
				
					//Call the function that executes the extract command
					exitStatus = calculation.runExtract();					
					
					//If the exit status of the function isn't 0, a problem occurred and we need to stop execution
					if (exitStatus != 0)
					{
						//Delete all files that may have been made during execution
						calculation.deleteOutput(0);
						
						//Display an error message
						JOptionPane.showMessageDialog(frame, "There was an error during the extract phase of execution.",
								"Filetype Error", JOptionPane.ERROR_MESSAGE);
						screenSwitch("menu");
						
						return null;
					}
					
					//Call the function that executes the nprofile command
					exitStatus = calculation.runNprofile();					
					
					//If the exit status of the function isn't 0, a problem occurred and we need to stop execution
					if (exitStatus != 0)
					{
						//Delete all files that may have been made during execution
						calculation.deleteOutput(0);
	
						//Display an error message
						JOptionPane.showMessageDialog(frame, "There was an error during the nprofile phase of execution.",
								"Filetype Error", JOptionPane.ERROR_MESSAGE);
						screenSwitch("menu");
						
						return null;
					}
					
					//Call the function that executes the acgt_gamma command
					calculation.runAcgt();
					
					//If the exit status of the function isn't 0, a problem occurred and we need to stop execution
					if (exitStatus != 0)
					{
						//Delete all files that may have been made during execution
						calculation.deleteOutput(0);
						
						//Display an error message
						JOptionPane.showMessageDialog(frame, "There was an error during the prediction phase of execution.",
								"Filetype Error", JOptionPane.ERROR_MESSAGE);
						screenSwitch("menu");
						
						return null;
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
							JOptionPane.showMessageDialog(frame, "There was an error during the plotting phase of execution.",
									"Filetype Error", JOptionPane.ERROR_MESSAGE);
							
							screenSwitch("menu");
							return null;
						}
						
						updateStatus("Generating page " + page_num + " out of " + allplotsResults[1]);
					}

					
					//Combine the results PS files into a single file
					calculation.combineOutput();

					//Delete the files used to make the combined PS file
					//Relocate the used files to an output directory
					calculation.cleanOutput();

					//statusLabel.setIcon(new ImageIcon(Toolkit.getDefaultToolkit().getImage("externals/img/complete.png").getScaledInstance(25, 25, Image.SCALE_DEFAULT)));
					statusLabel.setIcon(new ImageIcon(new ImageIcon(NpactGUI.class.getResource("img/complete.png")).getImage().getScaledInstance(25, 25, Image.SCALE_DEFAULT)));
					updateStatus("Analysis complete");
					directoryLabel1.setText("Results are located under folder: " + 
							"\"" + calculation.getResultFolder().toString() + "\"");
					directoryPanel.setVisible(true);
					returnPanel.setVisible(true);
					
					return null;
				}
			};
			task.execute();
		}
		//Handle "Show Plot" button press
		else if (event.getSource() == resultButton)
		{
			try
			{
				Desktop.getDesktop().open(new File(calculation.getNewDirectory()));
			}
			catch(Exception e)
			{
				JOptionPane.showMessageDialog(frame, "There was an error opening the results file",
						"File Open Error", JOptionPane.ERROR_MESSAGE);
			}
		}
		//Handle "New Analysis" button press
		else if (event.getSource() == returnButton)
		{
			screenSwitch("menu");
		}
	}
	
	//Function used to enable the configuration options for analysis
	private void enableConfig()
	{
		firstpgTextField.setText(calculation.getDescription());
		firstpgTextField.setEditable(true);
		followpgTextField.setText(calculation.getFollowingTitle());
		followpgTextField.setEditable(true);
		lengthTextField.setText(calculation.getLength());
		lengthTextField.setEditable(true);
		
		String nucleotides = calculation.getNucleotides();
		for (int i = 0; i < nucleotides.length(); i++)
		{
			if (nucleotides.charAt(i) == 'a')
			{
				aCBox.setSelected(true);
			}
			else if (nucleotides.charAt(i) == 'c')
			{
				cCBox.setSelected(true);
			}
			else if (nucleotides.charAt(i) == 'g')
			{
				gCBox.setSelected(true);
			}
			else if (nucleotides.charAt(i) == 't')
			{
				tCBox.setSelected(true);
			}
		}
		aCBox.setEnabled(true);
		cCBox.setEnabled(true);
		gCBox.setEnabled(true);
		tCBox.setEnabled(true);
		skipCBox.setEnabled(true);
		significanceBox.setEnabled(true);
		significanceBox.setSelectedIndex(1);
		sBaseTextField.setText("0");
		sBaseTextField.setEditable(true);
		eBaseTextField.setText(calculation.getLength());
		eBaseTextField.setEditable(true);
		startButton.setEnabled(true);
	}
	
	//Function used to disable the analysis configuration options
	private void disableConfig()
	{
		firstpgTextField.setText("");
		firstpgTextField.setEditable(false);
		followpgTextField.setText("");
		followpgTextField.setEditable(false);
		lengthTextField.setText("");
		lengthTextField.setEditable(false);
		aCBox.setSelected(false);
		aCBox.setEnabled(false);
		cCBox.setSelected(false);
		cCBox.setEnabled(false);
		gCBox.setSelected(false);
		gCBox.setEnabled(false);
		tCBox.setSelected(false);
		tCBox.setEnabled(false);
		significanceBox.setEnabled(false);
		significanceBox.setSelectedIndex(1);
		sBaseTextField.setText("");
		sBaseTextField.setEditable(false);
		eBaseTextField.setText("");
		eBaseTextField.setEditable(false);
		startButton.setEnabled(false);
	}
	
	//Call to update the text in the status label during program execution
	private void updateStatus(String text)
	{
		statusLabel.setText(text);
	}
	
	//Call to update the GUI to different screens
	private void screenSwitch(String text)
	{
		//Switch the GUI over the 'analysis in progress' information panels
		if (text.equals("working"))
		{
			containerPanel.removeAll();
			containerPanel.setLayout(new BoxLayout(containerPanel, BoxLayout.Y_AXIS));
			containerPanel.add(Box.createVerticalStrut(30));
			statusLabel.setIcon(new ImageIcon(new ImageIcon(NpactGUI.class.getResource("img/spinner.gif")).getImage().getScaledInstance(25, 25, Image.SCALE_DEFAULT)));
			statusLabel.setText("Analysis is being performed");
			
			directoryPanel.setVisible(false);
			returnPanel.setVisible(false);
			
			containerPanel.add(statusPanel);
			containerPanel.add(directoryPanel);
			containerPanel.add(Box.createVerticalStrut(30));
			containerPanel.add(returnPanel);
			containerPanel.add(Box.createVerticalStrut(220));
			containerPanel.revalidate();
			containerPanel.repaint();
		}
		//Switch the GUI over to the 'analysis configuration' information panels
		else if (text.equals("menu"))
		{
			containerPanel.removeAll();
			containerPanel.setLayout(new FlowLayout(FlowLayout.LEFT));
			containerPanel.add(boxPanel);
			
			disableConfig();
			
			calculation.resetGBKFile();
			
			fileTextField.setText("");
			containerPanel.revalidate();
			containerPanel.repaint();
		}
	}
}
