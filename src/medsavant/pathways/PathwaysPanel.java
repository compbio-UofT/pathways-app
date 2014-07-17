package medsavant.pathways;

import com.jidesoft.grid.SortableTable;
import com.jidesoft.pane.CollapsiblePane;
import com.jidesoft.swing.ButtonStyle;
import com.jidesoft.swing.JideButton;
import jannovar.common.VariantType;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Desktop;
import java.awt.Dimension;
import java.awt.Image;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.net.URI;
import java.net.URL;
import java.net.URLEncoder;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Properties;
import java.util.Set;
import javax.swing.*;
import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.table.DefaultTableColumnModel;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumn;
import net.miginfocom.swing.MigLayout;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.ut.biolab.medsavant.MedSavantClient;
import org.ut.biolab.medsavant.client.project.ProjectController;
import org.ut.biolab.medsavant.client.util.ClientMiscUtils;
import org.ut.biolab.medsavant.client.util.MedSavantWorker;
import org.ut.biolab.medsavant.client.view.component.ProgressWheel;
import org.ut.biolab.medsavant.client.view.component.SearchableTablePanel;
import org.ut.biolab.medsavant.client.view.dialog.IndividualSelector;
import org.ut.biolab.medsavant.client.view.genetics.variantinfo.ClinvarSubInspector;
import org.ut.biolab.medsavant.client.view.genetics.variantinfo.HGMDSubInspector;
import org.ut.biolab.medsavant.client.view.genetics.variantinfo.SimpleVariant;
import org.ut.biolab.medsavant.client.view.util.DialogUtils;
import org.ut.biolab.medsavant.client.view.util.ViewUtil;
import org.ut.biolab.medsavant.shared.appdevapi.AppColors;
import org.ut.biolab.medsavant.shared.format.AnnotationFormat;
import org.ut.biolab.medsavant.shared.format.BasicVariantColumns;
import org.ut.biolab.medsavant.shared.format.CustomField;

/**
 * Default panel for Pathways App.
 * @author rwong
 */
public class PathwaysPanel {
	
	private static final int SIDE_PANE_WIDTH= 380;
	
	private static Log log= LogFactory.getLog(MedSavantClient.class);
	private List<String> mutationFilterList= new LinkedList<String>();
	
	private final int PANE_WIDTH= 380;
	private final int PANE_WIDTH_OFFSET= 20;
        final String CACHEFOLDER = System.getProperty("user.home")+"/.medsavant/plugins/cache/";
        final String OUTPUTDIR = "pathways_plugin/";
        final String GPMLFOLDER = "gpml_files/";
        final String HTMLFOLDER = "html_files/";
	private String[] mutationArray;
	private JPanel appView= new JPanel();
	private JPanel optionsPanel= new JPanel();
	private JScrollPane scrollPane= new JScrollPane();
	private JPanel resultsPanel= new JPanel();
	private IndividualSelector patientSelector= new IndividualSelector(true);
	private JButton choosePatientButton;
        private JButton performAnalysisButton;
        private JButton openPathwayButton;
        private JLabel pathwayListLabel;
	private String currentHospitalID;
	private String currentDNAID;
        private JComboBox pathwayList;
        private String[][] pathwayInfo;
        private JTabbedPane tabbedPane;
        private JScrollPane tablePane;
        private JPanel tablePanel;
        private JScrollPane pngScrollPane;
        private JLabel pathwayImageLabel;
        private JFrame tableFrame;
        private SearchableTablePanel table;
        private Properties properties;
        private JScrollPane mutationPanelScroll;
        final private String PROPERTIESDIR = "properties/";
        final private String PROPERTIESFILENAME = "properties.xml";
        final private String DEFAULTPROPERTIESDIR = "medsavant/pathways/properties/default_properties.xml";
	
        final private Color LIGHTPURPLE = new Color(242, 202, 252);
        
        private CollapsiblePane mutationFilterGUI;
	final PathwayAnalysis pathwayAnalysisObject = new PathwayAnalysis();
	private static final List<String> JANNOVAR_MUTATIONS= Arrays.asList(
		VariantType.MISSENSE.toString(), VariantType.SYNONYMOUS.toString(),
		VariantType.FS_DELETION.toString(), VariantType.FS_INSERTION.toString(),
		VariantType.FS_SUBSTITUTION.toString(), VariantType.FS_DUPLICATION.toString(),
		VariantType.NON_FS_DELETION .toString(), VariantType.NON_FS_INSERTION.toString(),
		VariantType.NON_FS_SUBSTITUTION.toString(), VariantType.NON_FS_DUPLICATION.toString(),
		VariantType.SPLICING.toString(), VariantType.STOPGAIN.toString(),
		VariantType.STOPLOSS.toString(), 
		VariantType.START_LOSS.toString(), VariantType.UTR3.toString(),
		VariantType.UTR5.toString(), VariantType.INTRONIC.toString(),
		VariantType.UPSTREAM.toString(), VariantType.DOWNSTREAM.toString(),
		VariantType.INTERGENIC.toString(), VariantType.ncRNA_EXONIC.toString(),
		VariantType.ncRNA_INTRONIC.toString(), VariantType.ncRNA_SPLICING.toString(),
		VariantType.ERROR.toString()
	);
	private static final String[] DEFAULT_MUTATIONS= new String[] {
		VariantType.MISSENSE.toString(), VariantType.FS_DELETION.toString(),
		VariantType.FS_INSERTION.toString(), VariantType.FS_SUBSTITUTION.toString(),
		VariantType.FS_DUPLICATION.toString(), VariantType.NON_FS_DELETION .toString(),
		VariantType.NON_FS_INSERTION.toString(), VariantType.NON_FS_SUBSTITUTION.toString(),
		VariantType.NON_FS_DUPLICATION.toString(), VariantType.SPLICING.toString(),
		VariantType.STOPGAIN.toString(), VariantType.START_LOSS.toString()
	};
	private JLabel status;
	private ProgressWheel statusWheel;
	/**
	 * Create a new PathwaysPanel.
	 */
	public PathwaysPanel() {
		initView();
                
	}
	
	
	/**
	 * Return the main JPanel of this PathwaysPanel
	 * @return the main view JPanel
	 */
	public JPanel getView() {
		return appView;
	}
	
	
	/**
	 * Default initial view.
	 */
	private void initView() {
		// Create the options view
		optionsPanel.setLayout(new MigLayout("fillx"));
		optionsPanel.setMinimumSize(new Dimension(SIDE_PANE_WIDTH, 1));
		optionsPanel.setPreferredSize(new Dimension(SIDE_PANE_WIDTH, optionsPanel.getMaximumSize().height));
		optionsPanel.setBackground(LIGHTPURPLE);
		optionsPanel.setBorder(BorderFactory.createMatteBorder(0, 0, 0, 1, Color.LIGHT_GRAY));
		
		
		choosePatientButton= new JButton("Choose patient");
		choosePatientButton.addActionListener(choosePatientAction());
		
		optionsPanel.add(choosePatientButton, "alignx center, wrap");
		
                
                mutationFilterGUI = this.mutationCheckboxPanel();
                
                System.out.println("mutation width: "+mutationFilterGUI.getContentPaneWidth()+" height: "+mutationFilterGUI.getContentPaneHeight());
                System.out.println("Panel width: "+PANE_WIDTH+" offset: "+PANE_WIDTH_OFFSET);

                //mutationPanelScroll = new JScrollPane();
                //mutationPanelScroll.setViewportView(mutationFilterGUI);
                optionsPanel.add(mutationFilterGUI, "alignx center, wrap");
		
                optionsPanel.add(new JLabel(" "), "alignx center, wrap");
                
                performAnalysisButton = new JButton("Perform Hypergeometric Analysis");
                performAnalysisButton.addActionListener(analysisAction());
                performAnalysisButton.setVisible(false);
                optionsPanel.add(performAnalysisButton,"alignx center, wrap");
                
                pathwayListLabel = new JLabel("");
                optionsPanel.add(pathwayListLabel, "alignx center, wrap");
                
                status= new JLabel();
		statusWheel= new ProgressWheel();
		statusWheel.setIndeterminate(true);
		// hide for now
		status.setVisible(false);
		statusWheel.setVisible(false);
                
                pathwayList = new JComboBox(new String[0]);
                pathwayList.setPrototypeDisplayValue("this will be the width of the combo box");
                pathwayList.setVisible(false);
                optionsPanel.add(pathwayList,"alignx center, wrap");

                optionsPanel.add(status, "alignx center, wrap");
                optionsPanel.add(statusWheel, "alignx center, wrap");
                
                openPathwayButton = new JButton("Open selected pathway in browser");
                openPathwayButton.addActionListener(openPathwayAction());
                openPathwayButton.setVisible(false);
                optionsPanel.add(openPathwayButton, "alignx center, wrap");
                
		// Create the results view
		resultsPanel.setLayout(new MigLayout("fill"));
              
		resultsPanel.setBackground(Color.WHITE);
		
                resultPanelChoosePatientPrompt();
                
		// Add all components to our main view
                
		appView.setLayout(new MigLayout("insets 0 2 0 0, fillx, filly, wrap 1"));
                appView.setBackground(Color.LIGHT_GRAY);
		appView.add(optionsPanel,"dock west");
		appView.add(resultsPanel, "align center, growy, growx, push");//"growx, growy");
		// set the preferred size once the component is displayed.
		appView.addComponentListener(new ComponentListener()
			{				
				@Override
				public void componentShown(ComponentEvent e) {
					Dimension d= appView.getSize();
					scrollPane.setPreferredSize(new Dimension(d.width - SIDE_PANE_WIDTH, d.height));
					scrollPane.setMinimumSize(new Dimension(d.width - SIDE_PANE_WIDTH, d.height));
					scrollPane.setMaximumSize(new Dimension(d.width - SIDE_PANE_WIDTH, d.height));
					appView.updateUI();
				}
				
				@Override
				public void componentResized(ComponentEvent e) {
					componentShown(e);
				}
				
				@Override public void componentHidden(ComponentEvent e) {}
				@Override public void componentMoved(ComponentEvent e) {}
			}
		);
	}
	
	/**
	 * Action to perform when choose patient button is clicked.
	 * @return the ActionListener for this button
	 */
	private ActionListener choosePatientAction() {
		// create an anonymous class
		ActionListener outputAL= new ActionListener()
			{
				@Override
				public void actionPerformed(ActionEvent ae) {
					/* Show the patient selector window and get the patient selected
					 * by user. */
					patientSelector.setVisible(true);
					Set<String> selectedIndividuals= patientSelector.getHospitalIDsOfSelectedIndividuals();

					/* Once the user has made a patient hospital ID selection, get 
					 * the DNA ID so we can retrieve the patient's variants. */
					if (patientSelector.hasMadeSelection()) {
						currentHospitalID= patientSelector.getHospitalIDsOfSelectedIndividuals().iterator().next();
						String newDNAID= patientSelector.getDNAIDsOfSelectedIndividuals().iterator().next();

						if (newDNAID != null) {
							currentDNAID= newDNAID;
							choosePatientButton.setText(currentHospitalID);
                                                        performAnalysisButton.setVisible(true);
                                                        optionsPanel.revalidate();
                                                        
						} else { // can't find this individual's DNA ID - may be a DB error
							errorDialog("Can't find a DNA ID for " + currentHospitalID);
						}
					}

                                        resultPanelStartAnalysisPrompt();
				}
			};
		
		return outputAL;
	}
        
        /**
	 * Update the variantPane with the set of variants.
	 */
	private void updateResultsTable() {	
		table = pathwayAnalysisObject.getTableOutput();
                tablePane.setViewportView(table);
		
		table.scrollSafeSelectAction(new Runnable() {
            @Override
            public void run() {
				
                if (table.getTable().getSelectedRow() != -1) {
                    SortableTable st= PathwaysPanel.this.table.getTable();
                    int selectedIndex= st.getSelectedRow();
                    
                    String pathwayTitle = (String) st.getModel().getValueAt(selectedIndex, TestedPathway.PATHWAYNAMEINDEX);
                    
                    ImageIcon pngImage = pathwayAnalysisObject.getPathwayImage(pathwayTitle);
                    if (pngImage==null) {
                        pathwayImageLabel.setText("No image available for this pathway.");
                    }
                    else {
                        pathwayImageLabel.setIcon(pngImage);
                    }
                    saveNewProperties();
               }
            }
        });
	}
        
        private void resultPanelChoosePatientPrompt() {
            resultsPanel.removeAll();
            resultsPanel.add(new JLabel("Please choose patient for analysis"),"align center");
            resultsPanel.revalidate();
            resultsPanel.repaint();
        }
        
        private void resultPanelStartAnalysisPrompt() {
            resultsPanel.removeAll();
            resultsPanel.add(new JLabel("Please run analysis to see results"), "align center");
            resultsPanel.revalidate();
            resultsPanel.repaint();
            pathwayList.setVisible(false);
            status.setVisible(false);
            statusWheel.setVisible(false);
            openPathwayButton.setVisible(false);
        }
        
        private void resultsPanelAddTabbedPane() {
            resultsPanel.removeAll();
            tabbedPane = ViewUtil.getMSTabedPane();
            tabbedPane.setBackground(Color.CYAN);
            tablePane = new JScrollPane();
            tablePanel = new JPanel();
            tablePanel.setBackground(Color.PINK);
            tablePanel.setLayout(new BorderLayout());
            tablePanel.add(tablePane, BorderLayout.CENTER);
            tabbedPane.addTab("Results",tablePanel);
            pathwayImageLabel = new JLabel();
            pngScrollPane = new JScrollPane(pathwayImageLabel);
            tabbedPane.addTab("Pathway Visualization",pngScrollPane);
            tabbedPane.setToolTipTextAt(0,"P-values for the enrichment of each pathway");
            tabbedPane.setToolTipTextAt(1,"Image of selected pathway");
            resultsPanel.add(tabbedPane, "alignx center, growx, growy");
            resultsPanel.revalidate();
            
            resultsPanel.repaint();
        }
        
        /**
	 * Action to perform when choose patient button is clicked.
	 * @return the ActionListener for this button
	 */
	private ActionListener openPathwayAction() {
		// create an anonymous class
		ActionListener outputAL= new ActionListener()
			{
                            @Override
                            public void actionPerformed(ActionEvent ae) {
                                //generate HTML files
                                
                                /* Background task. */
                                MedSavantWorker generateVisualizationThread= new MedSavantWorker<Object>(PathwaysPanel.class.getCanonicalName()) {			
                                        @Override
                                        protected Object doInBackground() {

                                                /* Create and perform a new analysis. Uses a CountDownLatch to 
                                                 * ensure that currentPGXAnalysis is initilized before I can
                                                 * do anything with it (for example, cancel it). */
                                                try {
                                                        status.setText("Generating pathway visualization...");
                                                        status.setVisible(true);
                                                        System.out.println(pathwayInfo[2][0]+" and "+pathwayInfo[2][1]);
                                                        pathwayAnalysisObject.generateGPMLFile(pathwayInfo[0][pathwayList.getSelectedIndex()], pathwayInfo[2][pathwayList.getSelectedIndex()]);
                                                        pathwayAnalysisObject.generateHTMLFile(pathwayInfo[2][pathwayList.getSelectedIndex()]);
                                                        //cancelLatch.countDown();
                                                } catch (Exception e) {
                                                        errorDialog(e.getMessage());
                                                        e.printStackTrace();
                                                }

                                                return null;
                                        }

                                        @Override
                                        protected void showSuccess(Object t) {
                                            //add JComboBox list of pathways for selection
                                            status.setText("Finished generating pathway visualization");
                                            status.setVisible(true);
                                            statusWheel.setVisible(false);
                                            //open link
                                            
                                            String htmlPath = CACHEFOLDER+OUTPUTDIR+HTMLFOLDER;
                                            String pathStart = "file://";
                                            String[] test = pathwayInfo[1];
                                            System.out.println("test: "+test[0]+" and "+test[1]);
                                            System.out.println("attempting to open link at: "+pathStart+htmlPath+pathwayInfo[1][pathwayList.getSelectedIndex()]);
                                            openLink(pathStart+htmlPath+pathwayInfo[1][pathwayList.getSelectedIndex()]);
                                            System.out.println("link open attempted");
                                        }
                                };

                                // Execute thread
                                generateVisualizationThread.execute();  
                                
                                
                            }
			};
		
		return outputAL;
	}
        
        
    /**
     * Open link in browser.
     * @param htmlFileLink String link file path
     */
    private static void openLink(String htmlFileLink) {
        try {
            //URI uri = new URI(htmlFileLink);
            File file = new File(htmlFileLink);
            
            if (Desktop.isDesktopSupported()) {
              try {
                String operatingSystem = System.getProperty("os.name");
                if (operatingSystem.contains("Mac")) {
                    openUrlInBrowser(htmlFileLink);
                }
                else {
                    java.awt.Desktop.getDesktop().browse(file.toURI());
                    System.out.println("link open attempt complete");
                }
                
              } catch (IOException e) { e.printStackTrace(); }
            } else { System.out.println("desktop not supported - cannot open link in browsers"); }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private static void openUrlInBrowser(String url)
    {
        Runtime runtime = Runtime.getRuntime();
        String[] args = { "osascript", "-e", "open location \"" + url + "\"" };
        try
        {
            Process process = runtime.exec(args);
        }
        catch (IOException e)
        {
            e.printStackTrace();
        }
    }
      
        
        /**
	 * Action to perform when choose patient button is clicked.
	 * @return the ActionListener for this button
	 */
	private ActionListener analysisAction() {
		// create an anonymous class
		ActionListener outputAL= new ActionListener()
			{
				@Override
				public void actionPerformed(ActionEvent ae) {
				//show status message and progress wheel
                                status.setText("Analyzing pathway enrichment...");
                                status.setVisible(true);
                                statusWheel.setVisible(true);
                                    
                /* Background task. */
		MedSavantWorker pathwayAnalysisThread= new MedSavantWorker<Object>(PathwaysPanel.class.getCanonicalName()) {
                        @Override
			protected Object doInBackground() {
				
                                choosePatientButton.setEnabled(false);
				/* Create and perform a new analysis. Uses a CountDownLatch to 
				 * ensure that currentPGXAnalysis is initilized before I can
				 * do anything with it (for example, cancel it). */
				try {
                                        System.out.println("CURRENT DNAID IS: "+currentDNAID);
					pathwayInfo = pathwayAnalysisObject.hypergeometricWithWikiPathways(currentDNAID, mutationFilterList);
					//cancelLatch.countDown();
				} catch (Exception e) {
					errorDialog(e.getMessage());
					e.printStackTrace();
				}

				return null;
			}

			@Override
			protected void showSuccess(Object t) {
                            resultsPanelAddTabbedPane();
                            choosePatientButton.setEnabled(true);
                            status.setText("Finished analyzing pathway enrichment");
                            status.setVisible(true);
                            statusWheel.setVisible(false);
                            //make combo box with pathway titles for user to select from
                            //DELETE LATER - DO EVERYTHING FROM TABLE
                            addPathwayTitlesToSelectionBox(pathwayInfo[0], pathwayInfo[3]);
                            pathwayList.setVisible(true);
                            //add OK button
                            openPathwayButton.setVisible(true);
                            optionsPanel.getRootPane().revalidate();
                            updateResultsTable();
                            setTableColumnWidths();
                            resultsPanel.revalidate();
                            
			}
		};
		
		// Execute thread
		pathwayAnalysisThread.execute();  
                                    

					/* Prevent further patient selection while an analysis thread is
					 * running. */
					//choosePatientButton.setEnabled(false);

					/* Perform a analysis. */
					// CALL ANALYSIS METHOD IN A NEW THREAD
				}
			};
		
		return outputAL;
	}
        
        /**
	 * Fill combobox with pathway titles for selection
	 */
        private void addPathwayTitlesToSelectionBox(String[] pathwayTitles, String[] pValues) {
            for (int i = 0; i < pathwayTitles.length; i++) {
                if (pValues[i].equals("1")) {
                    pathwayList.addItem(pathwayTitles[i]);
                }
                else {
                    pathwayList.addItem(pValues[i]+" "+pathwayTitles[i]);
                }
                
            }
        }
        
	/**
	 * Create an error dialog and output the error to the log.
	 * @param errorMessage the error message to display.
	 */
	private void errorDialog(String errorMessage) {
		DialogUtils.displayError("Oops!", errorMessage);
		log.error("[" + this.getClass().getSimpleName() + "]: " + errorMessage);
	}
	
        /**
	 * Create a CollapsiblePane of a checkbox panel of mutations
	 * @return A mutation checkbox CollapsiblePane
	 */
	private CollapsiblePane mutationCheckboxPanel() {
		CollapsiblePane collapsibleMutation= new CollapsiblePane("Mutations");
		collapsibleMutation.setLayout(new MigLayout("gapy 0px"));
		
                JPanel mutationCheckList = new JPanel();
                mutationCheckList.setLayout(new MigLayout("gapy 0px"));
		for (String jm : JANNOVAR_MUTATIONS) {
			final JCheckBox currentCheckBox= new JCheckBox(jm);
			
			// Allow checkboxes to register themselves as checked or unchecked upon being clicked
			currentCheckBox.addActionListener(
				new ActionListener() {
					@Override
					public void actionPerformed(ActionEvent e) {
						if (currentCheckBox.isSelected()) {
							mutationFilterList.add(currentCheckBox.getText());
						} else {
							mutationFilterList.remove(currentCheckBox.getText());
						}
					}
				}
			);
			
			// Set the defaults
			if (Arrays.asList(DEFAULT_MUTATIONS).contains(jm)) {
				currentCheckBox.setSelected(true);
				mutationFilterList.add(currentCheckBox.getText());
			}
			
			mutationCheckList.add(currentCheckBox, "wrap");
		}
		mutationCheckList.setBackground(LIGHTPURPLE);
		collapsibleMutation.setStyle(CollapsiblePane.PLAIN_STYLE);
		collapsibleMutation.setFocusPainted(false);
		collapsibleMutation.collapse(true);	
		
		collapsibleMutation.setMinimumSize(new Dimension(PANE_WIDTH - PANE_WIDTH_OFFSET, 0));
		collapsibleMutation.setContentPaneHeight(450);
                
                //add scroll
                JScrollPane scrollPane = new JScrollPane();
                scrollPane.setViewportView(mutationCheckList);
                collapsibleMutation.setContentPane(scrollPane);
		return collapsibleMutation;
	}
        
    private void getProperties() {
        //check if properties dir exists
        properties = new Properties();
        File propDir = new File(CACHEFOLDER+OUTPUTDIR+PROPERTIESDIR);
        if (!propDir.exists()) {
            propDir.mkdir();
        }
        else if ((new File(CACHEFOLDER+OUTPUTDIR+PROPERTIESDIR+PROPERTIESFILENAME)).exists()) {
            try {
                properties.loadFromXML(new FileInputStream(new File(CACHEFOLDER+OUTPUTDIR+PROPERTIESDIR+PROPERTIESFILENAME)));
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }
        else {
            try {
                properties.loadFromXML(this.getClass().getClassLoader().getResourceAsStream(DEFAULTPROPERTIESDIR));
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }
    }
    
    private void setTableColumnWidths() {
        this.getProperties();
        SortableTable sortableTable = table.getTable();
        sortableTable.getColumnModel().getColumn(TestedPathway.PVALUEINDEX).setPreferredWidth(Integer.parseInt(properties.getProperty("pValueColumnWidth")));
        sortableTable.getColumnModel().getColumn(TestedPathway.PATHWAYNAMEINDEX).setPreferredWidth(Integer.parseInt(properties.getProperty("pathwayNameColumnWidth")));
        sortableTable.getColumnModel().getColumn(TestedPathway.GENESINDEX).setPreferredWidth(Integer.parseInt(properties.getProperty("genesColumnWidth")));
        
    }
    
    private void saveNewProperties() {
        SortableTable sortableTable = table.getTable();
        properties.setProperty("pValueColumnWidth",sortableTable.getColumnModel().getColumn(TestedPathway.PVALUEINDEX).getWidth()+"");
        properties.setProperty("pathwayNameColumnWidth",sortableTable.getColumnModel().getColumn(TestedPathway.PATHWAYNAMEINDEX).getWidth()+"");
        properties.setProperty("genesColumnWidth",sortableTable.getColumnModel().getColumn(TestedPathway.GENESINDEX).getWidth()+"");
        File newProperties = new File(CACHEFOLDER+OUTPUTDIR+PROPERTIESDIR+PROPERTIESFILENAME);
        
        try {
            if (!newProperties.exists()) {
                newProperties.createNewFile();
            }
            properties.storeToXML(new FileOutputStream(newProperties),"custom column widths from the last time the user generated a pathway visualization");
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    
    /**
	 * Create a button that opens a URL in a web browser when clicked.
	 * @param buttonText Button text
	 * @param baseURL URL linked from the button
	 * @param appendToURL Append text to the URL
	 * @param doEncode encode the text using the UTF-8
	 * @return a JideButton that opens the URL in a web browser
	 */
	private JideButton getURLButton(String buttonText, final String baseURL, 
		final String appendToURL, final boolean doEncode) {
		
		final String URL_CHARSET = "UTF-8";
		
		JideButton urlButton= new JideButton(buttonText);
		urlButton.setButtonStyle(ButtonStyle.HYPERLINK_STYLE);
		urlButton.setForeground(Color.BLUE);
		urlButton.setToolTipText("Lookup " + buttonText + " on the web");
		urlButton.addActionListener(new ActionListener() 
		{
			@Override
			public void actionPerformed(ActionEvent ae) {
				try {
					URL url;
					if (doEncode)
						url = new URL(baseURL + URLEncoder.encode(appendToURL, URL_CHARSET));
					else
						url = new URL(baseURL + appendToURL);
					
					java.awt.Desktop.getDesktop().browse(url.toURI());
				} catch (Exception ex) {
					ClientMiscUtils.reportError("Problem launching website: %s", ex);
				}
			}
		});
		
		return urlButton;
	}
	
	
	
}








