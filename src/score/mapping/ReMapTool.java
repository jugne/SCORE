package score.mapping;

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.Font;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.swing.BoxLayout;
import javax.swing.GroupLayout;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;
import javax.swing.WindowConstants;
import javax.swing.border.EtchedBorder;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.xml.sax.SAXException;

import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import feast.fileio.logfileiterator.LogFileRealParameter;
import feast.fileio.logfileiterator.TraceLogFileState;
import feast.function.Concatenate;
import score.logger.TypedNetworkStatsLogger;
import score.networkAnnotator.ExtendedNetworkBuilder;
import score.utils.LogFileIterator;
import score.utils.NetworkLogFileState;

/**
 * @author Ugne Stolz <ugne.stolz@protonmail.com>
 * @date 23 Apr 2020
 */
public class ReMapTool {


    private static class ReMappOptions {
		File xmlFile;
		File logFile;
		File networkFile;
		File outFile;
		Integer dimension;
		Boolean rejection;

        @Override
        public String toString() {
            return "Active options:\n" +
					"SCoRe XML file: " + xmlFile + "\n" +
					"Network file: " + networkFile + "\n" +
					"Log file: " + logFile + "\n" +
					"Output file: " + outFile +"\n" +
					"Parameters dimension: " + dimension +"\n"+
					"Rejection mapping: " + rejection;
        }
    }

	public ReMapTool(ReMappOptions options) throws Exception, ParserConfigurationException, SAXException, IOException {

        // Display options:
        System.out.println(options + "\n");

		BufferedReader reader = new BufferedReader(new FileReader(options.logFile));
		String header = reader.readLine();
		while (header.charAt(0) == '#'){
			header = reader.readLine();
		}
		reader.close();

		String[] headerSplit = header.split("\\t");
		Pattern migPattern = Pattern.compile("b_migrationRate.*");
		Pattern reaPattern = Pattern.compile("reassortmentRate.*");
		Pattern NePattern = Pattern.compile("popSize.*");

		// get header strings for all logger types
		List<String> migrationFromTo = new ArrayList<>();
		List<String> reassortment = new ArrayList<>();
		List<String> Ne = new ArrayList<>();

		for (String param : headerSplit) {
			Matcher migMatcher = migPattern.matcher(param);
			Matcher reaMatcher = reaPattern.matcher(param);
			Matcher NeMatcher = NePattern.matcher(param);
			if (migMatcher.find())
				migrationFromTo.add(migMatcher.group(0));
			if (reaMatcher.find())
				reassortment.add(reaMatcher.group(0));
			if (NeMatcher.find())
				Ne.add(NeMatcher.group(0));
		}

		if (migrationFromTo.size() == 0)
			throw new Exception("Backward migration not found in log file.");
		if (reassortment.size() == 0)
			throw new Exception("Reassortment not found in log file.");
		if (Ne.size() == 0)
			throw new Exception("PopSize not found in log file.");

		if (Ne.size()>options.dimension){
			List<String> Ne_ = new ArrayList<>();
			for (int i=0; i< options.dimension; i++){
				Ne_.add(Ne.get(i));
			}
			Ne = Ne_;
		}

		if (reassortment.size()>options.dimension){
			List<String> reassortment_ = new ArrayList<>();
			for (int i=0; i< options.dimension; i++){
				reassortment_.add(Ne.get(i));
			}
			reassortment = reassortment_;
		}


		if (Ne.size() != reassortment.size())
			throw new Exception("PopSize and Reassortment rate have different dimensions.");
		if (migrationFromTo.size() != reassortment.size() * (reassortment.size() - 1))
			throw new Exception("Backwards migration rate dimension: " + migrationFromTo.size() +
					"\nExpected: " + reassortment.size() * (reassortment.size() - 1));

		DocumentBuilderFactory docFactory = DocumentBuilderFactory.newInstance();
		DocumentBuilder docBuilder = docFactory.newDocumentBuilder();
		Document scoreXML = docBuilder.parse(options.xmlFile);
        
		// retrieve the element 'link'
		Element run = (Element) scoreXML.getElementsByTagName("run").item(0);

		// get parent node of run which will be used later to create a new run node
		Element beast = (Element) run.getParentNode();

		// remove the old run node
		beast.removeChild(run);

		Element feastRun = scoreXML.createElement("run");
		beast.appendChild(feastRun);
		feastRun.setAttribute("spec", LogFileIterator.class.getCanonicalName());

		// trace log state
		Element traceLogState = scoreXML.createElement("logFileState");
		traceLogState.setAttribute("spec", TraceLogFileState.class.getCanonicalName());
		traceLogState.setAttribute("logFileName", options.logFile.getName());
		feastRun.appendChild(traceLogState);

		List<Element> reaLogFileEntries = new ArrayList<>();
		List<Element> reafieldParameters = new ArrayList<>();
		for (int i = 0; i < reassortment.size(); i++) {
			reaLogFileEntries.add(scoreXML.createElement("logFileEntry"));
			reaLogFileEntries.get(i).setAttribute("spec", LogFileRealParameter.class.getCanonicalName());
			reaLogFileEntries.get(i).setAttribute("fieldName", reassortment.get(i));

			reafieldParameters.add(scoreXML.createElement("fieldParameter"));
			reafieldParameters.get(i).setAttribute("id", reassortment.get(i));
			reafieldParameters.get(i).setAttribute("spec", RealParameter.class.getCanonicalName());
			reafieldParameters.get(i).setAttribute("value", "0.0");

			reaLogFileEntries.get(i).appendChild(reafieldParameters.get(i));
			traceLogState.appendChild(reaLogFileEntries.get(i));
		}

		List<Element> NeLogFileEntries = new ArrayList<>();
		List<Element> NefieldParameters = new ArrayList<>();
		for (int i = 0; i < reassortment.size(); i++) {
			NeLogFileEntries.add(scoreXML.createElement("logFileEntry"));
			NeLogFileEntries.get(i).setAttribute("spec", LogFileRealParameter.class.getCanonicalName());
			NeLogFileEntries.get(i).setAttribute("fieldName", Ne.get(i));

			NefieldParameters.add(scoreXML.createElement("fieldParameter"));
			NefieldParameters.get(i).setAttribute("id", Ne.get(i));
			NefieldParameters.get(i).setAttribute("spec", RealParameter.class.getCanonicalName());
			NefieldParameters.get(i).setAttribute("value", "0.0");

			NeLogFileEntries.get(i).appendChild(NefieldParameters.get(i));
			traceLogState.appendChild(NeLogFileEntries.get(i));
		}

		List<Element> migLogFileEntries = new ArrayList<>();
		List<Element> migfieldParameters = new ArrayList<>();
		for (int i = 0; i < reassortment.size(); i++) {
			migLogFileEntries.add(scoreXML.createElement("logFileEntry"));
			migLogFileEntries.get(i).setAttribute("spec", LogFileRealParameter.class.getCanonicalName());
			migLogFileEntries.get(i).setAttribute("fieldName", migrationFromTo.get(i));

			migfieldParameters.add(scoreXML.createElement("fieldParameter"));
			migfieldParameters.get(i).setAttribute("id", migrationFromTo.get(i));
			migfieldParameters.get(i).setAttribute("spec", RealParameter.class.getCanonicalName());
			migfieldParameters.get(i).setAttribute("value", "0.0");

			migLogFileEntries.get(i).appendChild(migfieldParameters.get(i));
			traceLogState.appendChild(migLogFileEntries.get(i));
		}

		// network log state
		Element networkLogState = scoreXML.createElement("logFileState");
		networkLogState.setAttribute("spec", NetworkLogFileState.class.getCanonicalName());
		networkLogState.setAttribute("logFileName", options.networkFile.getName());

		Element network = scoreXML.createElement("network");
		network.setAttribute("spec", ExtendedNetworkBuilder.class.getCanonicalName());
		network.setAttribute("id", "network");
		networkLogState.appendChild(network);

		Element taxonSet = scoreXML.createElement("taxonSet");
		taxonSet.setAttribute("idref", "taxonSet");
		networkLogState.appendChild(taxonSet);

		feastRun.appendChild(networkLogState);

		// loggers

		// network logger
		Element netLogger = scoreXML.createElement("logger");
		netLogger.setAttribute("spec", "Logger");
		netLogger.setAttribute("logEvery", "1");
		netLogger.setAttribute("mode", "tree");
		netLogger.setAttribute("fileName", "$(filebase).network.trees");

		Element netLog = scoreXML.createElement("log");
		netLog.setAttribute("spec", MappedNetwork.class.getCanonicalName());
		netLog.setAttribute("untypedNetwork", "@network");
		netLog.setAttribute("dimension", Integer.toString(migLogFileEntries.size()));
		if (options.rejection !=null)
			netLog.setAttribute("rejection", options.rejection.toString());
		netLog.setAttribute("mapOnInit", "false");
		netLog.setAttribute("id", "mappedNet");

		Element NePopSize = scoreXML.createElement("Ne");
		NePopSize.setAttribute("id", "popSize");
		NePopSize.setAttribute("spec", Concatenate.class.getCanonicalName());
		List<Element> NeArg = new ArrayList<>();
		for (int i = 0; i < NeLogFileEntries.size(); i++) {
			NeArg.add(scoreXML.createElement("arg"));
//			NeArg.get(i).setAttribute("spec", "RealParameter");
			NeArg.get(i).setAttribute("idref", Ne.get(i));
			NePopSize.appendChild(NeArg.get(i));
		}
		netLog.appendChild(NePopSize);

		Element backwardsMig = scoreXML.createElement("backwardsMigration");
		backwardsMig.setAttribute("id", "migrationRate");
		backwardsMig.setAttribute("spec", Concatenate.class.getCanonicalName());
		List<Element> backwardsMigArg = new ArrayList<>();
		for (int i = 0; i < migLogFileEntries.size(); i++) {
			backwardsMigArg.add(scoreXML.createElement("arg"));
//			backwardsMigArg.get(i).setAttribute("spec", "RealParameter");
			backwardsMigArg.get(i).setAttribute("idref", migrationFromTo.get(i));
			backwardsMig.appendChild(backwardsMigArg.get(i));
		}
		netLog.appendChild(backwardsMig);

		Element reassortmentRates = scoreXML.createElement("reassortmentRates");
		reassortmentRates.setAttribute("id", "reassortmentRate");
		reassortmentRates.setAttribute("spec", Concatenate.class.getCanonicalName());
		List<Element> reassortmentRatesArg = new ArrayList<>();
		for (int i = 0; i < reaLogFileEntries.size(); i++) {
			reassortmentRatesArg.add(scoreXML.createElement("arg"));
//			reassortmentRatesArg.get(i).setAttribute("spec", "RealParameter");
			reassortmentRatesArg.get(i).setAttribute("idref", reassortment.get(i));
			reassortmentRates.appendChild(reassortmentRatesArg.get(i));
		}
		netLog.appendChild(reassortmentRates);

		Element typeTrait = scoreXML.createElement("typeTrait");
		typeTrait.setAttribute("idref", "typeTrait");
		netLog.appendChild(typeTrait);
		netLogger.appendChild(netLog);
		feastRun.appendChild(netLogger);

		// stat logger
		Element mapStat = scoreXML.createElement("logger");
		mapStat.setAttribute("spec", "Logger");
		mapStat.setAttribute("logEvery", "1");
		mapStat.setAttribute("fileName", "$(filebase).stats.log");

		Element netStatLog = scoreXML.createElement("log");
		netStatLog.setAttribute("spec", TypedNetworkStatsLogger.class.getCanonicalName());
		netStatLog.setAttribute("network", "@mappedNet");

		mapStat.appendChild(netStatLog);
		feastRun.appendChild(mapStat);

		// create the xml file
		// transform the DOM Object to an XML File
		TransformerFactory transformerFactory = TransformerFactory.newInstance();
		Transformer transformer = transformerFactory.newTransformer();
		transformer.setOutputProperty(OutputKeys.INDENT, "yes");
		transformer.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "4");
		DOMSource domSource = new DOMSource(scoreXML);
		StreamResult streamResult = new StreamResult(options.outFile);

		transformer.transform(domSource, streamResult);

		System.out.println("Done creating re-mapp XML File");

    }

        
    /**
     * Use a GUI to retrieve ACGAnnotator options.
     *
     * @param options options object to populate using GUI
     * @return true if options successfully collected, false otherwise
     */
    private static boolean getOptionsGUI(ReMappOptions options) {

        boolean[] canceled = {false};

        JDialog dialog = new JDialog((JDialog)null, true);
        dialog.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        dialog.setLocationRelativeTo(null);
		dialog.setTitle("SCoRe Re-Mapper");

		JLabel xmlFileLabel = new JLabel("SCoRe XML file:");
		JLabel netFileLabel = new JLabel("Network file:");
		JLabel logFileLabel = new JLabel("Log file:");
        JLabel outFileLabel = new JLabel("Output file:");

		JTextField xmlFilename = new JTextField(20);
		xmlFilename.setEditable(false);
		JButton xmlFileButton = new JButton("Choose File");

		JTextField netFilename = new JTextField(20);
		netFilename.setEditable(false);
		JButton netFileButton = new JButton("Choose File");

		JTextField logFilename = new JTextField(20);
		logFilename.setEditable(false);
		JButton logFileButton = new JButton("Choose File");

        JTextField outFilename = new JTextField(20);
        outFilename.setEditable(false);
        JButton outFileButton = new JButton("Choose File");


        Container cp = dialog.getContentPane();
        BoxLayout boxLayout = new BoxLayout(cp, BoxLayout.PAGE_AXIS);
        cp.setLayout(boxLayout);

        JPanel mainPanel = new JPanel();

        GroupLayout layout = new GroupLayout(mainPanel);
        mainPanel.setLayout(layout);
        layout.setAutoCreateGaps(true);
        layout.setAutoCreateContainerGaps(true);

        layout.setHorizontalGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup()
						.addComponent(xmlFileLabel)
						.addComponent(netFileLabel)
						.addComponent(logFileLabel)
						.addComponent(outFileLabel))
				.addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
						.addComponent(xmlFilename)
						.addComponent(netFilename)
						.addComponent(logFilename)
						.addComponent(outFilename))
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
						.addComponent(xmlFileButton)
						.addComponent(netFileButton)
						.addComponent(logFileButton)
                        .addComponent(outFileButton))
                );

        layout.setVerticalGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup()
						.addComponent(xmlFileLabel)
						.addComponent(xmlFilename,
								GroupLayout.PREFERRED_SIZE,
								GroupLayout.DEFAULT_SIZE,
								GroupLayout.PREFERRED_SIZE)
						.addComponent(xmlFileButton))
				.addGroup(layout.createParallelGroup()
						.addComponent(netFileLabel)
						.addComponent(netFilename,
								GroupLayout.PREFERRED_SIZE,
								GroupLayout.DEFAULT_SIZE,
								GroupLayout.PREFERRED_SIZE)
						.addComponent(netFileButton))
				.addGroup(layout.createParallelGroup()
                        .addComponent(logFileLabel)
						.addComponent(logFilename,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE)
						.addComponent(logFileButton))
                .addGroup(layout.createParallelGroup()
                        .addComponent(outFileLabel)
                        .addComponent(outFilename,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE)
                        .addComponent(outFileButton))
                );

        mainPanel.setBorder(new EtchedBorder());
        cp.add(mainPanel);

        JPanel buttonPanel = new JPanel();

        JButton runButton = new JButton("Analyze");
		runButton.addActionListener((e) -> {
			dialog.setVisible(false);
		});
        runButton.setEnabled(false);
        buttonPanel.add(runButton);

        JButton cancelButton = new JButton("Quit");
        cancelButton.addActionListener((e) -> {
            dialog.setVisible(false);
            canceled[0] = true;
        });
        buttonPanel.add(cancelButton);

		JFileChooser xmlFileChooser = new JFileChooser();
		xmlFileButton.addActionListener(e -> {
			xmlFileChooser.setDialogTitle("Select SCoRe XML file");
			if (options.xmlFile == null)
				xmlFileChooser.setCurrentDirectory(new File(System.getProperty("user.dir")));
			int returnVal = xmlFileChooser.showOpenDialog(dialog);

			if (returnVal == JFileChooser.APPROVE_OPTION) {
				options.xmlFile = xmlFileChooser.getSelectedFile();
				xmlFilename.setText(xmlFileChooser.getSelectedFile().getName());
				runButton.setEnabled(true);
			}
		});

		JFileChooser netFileChooser = new JFileChooser();
		netFileButton.addActionListener(e -> {
			netFileChooser.setDialogTitle("Select SCoRe Network file to remap");
            if (options.networkFile == null)
				netFileChooser.setCurrentDirectory(options.xmlFile);
			else
				netFileChooser.setCurrentDirectory(new File(System.getProperty("user.dir")));
			int returnVal = netFileChooser.showOpenDialog(dialog);

			if (returnVal == JFileChooser.APPROVE_OPTION) {
				options.networkFile = netFileChooser.getSelectedFile();
				netFilename.setText(netFileChooser.getSelectedFile().getName());
				runButton.setEnabled(true);
			}
		});

		JFileChooser logFileChooser = new JFileChooser();
		logFileButton.addActionListener(e -> {
			logFileChooser.setDialogTitle("Select SCoRe log file");
			if (options.logFile == null)
				logFileChooser.setCurrentDirectory(options.xmlFile);
			else
				logFileChooser.setCurrentDirectory(new File(System.getProperty("user.dir")));
			int returnVal = logFileChooser.showOpenDialog(dialog);

            if (returnVal == JFileChooser.APPROVE_OPTION) {
				options.logFile = logFileChooser.getSelectedFile();
				logFilename.setText(logFileChooser.getSelectedFile().getName());
                runButton.setEnabled(true);
            }
        });

        JFileChooser outFileChooser = new JFileChooser();
        outFileButton.addActionListener(e -> {
			outFileChooser.setDialogTitle("Select output XML file name");
            if (options.xmlFile != null)
                outFileChooser.setCurrentDirectory(options.xmlFile);
            else
                outFileChooser.setCurrentDirectory(new File(System.getProperty("user.dir")));

//            outFileChooser.setSelectedFile(options.outFile);
            int returnVal = outFileChooser.showOpenDialog(dialog);

            if (returnVal == JFileChooser.APPROVE_OPTION) {
                options.outFile = outFileChooser.getSelectedFile();
                outFilename.setText(outFileChooser.getSelectedFile().getName());
				runButton.setEnabled(true);
            }
        });

        cp.add(buttonPanel);

        dialog.pack();
        dialog.setResizable(false);
        dialog.setVisible(true);

        return !canceled[0];
    }

    /**
     * Prepare JFrame to which ACGAnnotator output streams will be
     * directed.
     */
    private static void setupGUIOutput() {

        JFrame frame = new JFrame();
        frame.setTitle("Reassortment Event Locator");
        frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);

        JTextArea textArea = new JTextArea(25, 80);
        textArea.setFont(new Font("monospaced", Font.PLAIN, 12));
        textArea.setEditable(false);
        frame.getContentPane().add(new JScrollPane(textArea), BorderLayout.CENTER);

        JButton closeButton = new JButton("Close");
        closeButton.addActionListener(e -> System.exit(0));
        JPanel buttonPanel = new JPanel();
        buttonPanel.add(closeButton);
        frame.getContentPane().add(buttonPanel, BorderLayout.PAGE_END);

        // Redirect streams to output window:
        OutputStream out = new OutputStream() {
            @Override
            public void write(int b) throws IOException {
                SwingUtilities.invokeLater(() -> {
                    if ((char)b == '\r') {
                        int from = textArea.getText().lastIndexOf("\n") + 1;
                        int to = textArea.getText().length();
                        textArea.replaceRange(null, from, to);
                    } else
                        textArea.append(String.valueOf((char) b));
                });
            }
        };

        System.setOut(new PrintStream(out, true));
        System.setErr(new PrintStream(out, true));

        frame.pack();
        frame.setVisible(true);
    }

    public static String helpMessage =
            "TrunkReassortment - counts how many reassortment events happened on trunk and non-trunk nodes.\n"
                    + "\n"
                    + "Usage: appstore ACGAnnotator [-help | [options] logFile [outputFile]\n"
                    + "\n"
                    + "Option                   Description\n"
                    + "--------------------------------------------------------------\n"
                    + "-help                    Display usage info.\n"
					+ "-xml 					XML file of the original SCoRe analysis.\n"
					+ "-network				    network.trees file produced by ScoRe.\n"
					+ "-log     				log file produced by ScoRe.\n"
					+ "-out     				Name of the re-mapping XML file.\n"
					+ "                         Mapped network log file names will follow\n"
					+ "                         this file's naming scheme.\n"
					+ "\n";

    /**
     * Print usage info and exit.
     */
    public static void printUsageAndExit() {
        System.out.println(helpMessage);
        System.exit(0);
    }

    /**
     * Display error, print usage and exit with error.
     */
    public static void printUsageAndError(String errMsg) {
        System.err.println(errMsg);
        System.err.println(helpMessage);
        System.exit(1);
    }

    /**
     * Retrieve TrunkReassortment options from command line.
     *
     * @param args command line arguments
     * @param options object to populate with options
     */
	public static void getCLIOptions(String[] args, ReMappOptions options) {
        int i=0;
		while (args.length > i && args[i].startsWith("-")) {
            switch(args[i]) {
                case "-help":
                    printUsageAndExit();
                    break;

			case "-xml":
                    if (args.length<=i+1)
					printUsageAndError("-xml must be followed by a xml file path.");

                    try {
					options.xmlFile = new File(args[i + 1]);
                    } catch (NumberFormatException e) {
					printUsageAndError("Error parsing XML file.");
                    }

                    i += 1;
                    break;


			case "-network":
                    if (args.length<=i+1) {
					printUsageAndError("-network must be followed by a network file path.");
                    }

                    try {
					options.networkFile = new File(args[i + 1]);
                    } catch (NumberFormatException e) {
					printUsageAndError("Error parsing network file.");
				}

                    i += 1;
                    break;
                    
			case "-log":
				if (args.length <= i + 1) {
					printUsageAndError("-log must be followed by a log file path.");
				}

				try {
					options.logFile = new File(args[i + 1]);
				} catch (NumberFormatException e) {
					printUsageAndError("Error parsing log file.");
				}

				i += 1;
				break;

			case "-out":
                    if (args.length<=i+1) {
					printUsageAndError("-out must be followed by an output file path.");
                    }

				try {
					options.outFile = new File(args[i + 1]);
				} catch (NumberFormatException e) {
					printUsageAndError("Error parsing output file path.");
				}

                    i += 1;
                    break;
			case "-dim":
					if (args.length<=i+1) {
						printUsageAndError("-out must be followed by an output file path.");
					}

					try {
						options.dimension = new Integer(args[i + 1]);
					} catch (NumberFormatException e) {
						printUsageAndError("Error parsing output file path.");
					}

					i += 1;
					break;
			case "-reject":
					if (args.length<=i+1) {
						printUsageAndError("-out must be followed by an output file path.");
					}

					try {
						options.rejection = new Boolean(args[i + 1]);
					} catch (NumberFormatException e) {
						printUsageAndError("Error parsing output file path.");
					}

					i += 1;
					break;


                default:
                    printUsageAndError("Unrecognised command line option '" + args[i] + "'.");
            }

            i += 1;
        }
    }

    /**
     * Main method for ACGAnnotator.  Sets up GUI if needed then
     * uses the ACGAnnotator constructor to actually perform the analysis.
     *
     * @param args command line arguments
     */
    public static void main(String[] args) {
		ReMappOptions options = new ReMappOptions();

        if (args.length == 0) {
            // Retrieve options from GUI:

            try {
                UIManager.setLookAndFeel(UIManager.getCrossPlatformLookAndFeelClassName());
            } catch (ClassNotFoundException | InstantiationException | UnsupportedLookAndFeelException | IllegalAccessException e) {
                Log.warning.println("Error setting cross-platform look and feel.");
            }

            try {
                SwingUtilities.invokeAndWait(() -> {
                    if (!getOptionsGUI(options))
                        System.exit(0);

                    setupGUIOutput();
                });
            } catch (InterruptedException | InvocationTargetException e) {
                e.printStackTrace();
            }


        } else {
            getCLIOptions(args, options);
        }

        // Run ACGAnnotator
        try {
			new ReMapTool(options);

        } catch (Exception e) {
            if (args.length == 0) {
                JOptionPane.showMessageDialog(null, e.getMessage(),
                        "Error", JOptionPane.ERROR_MESSAGE);
            } else {
                System.err.println("Error: " + e.getMessage());
                e.printStackTrace();
                System.err.println();
                System.err.println(helpMessage);
            }

            System.exit(1);
        }
    }
}