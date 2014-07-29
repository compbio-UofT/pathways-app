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
import java.awt.Font;
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
import java.nio.file.FileSystem;
import java.nio.file.FileSystems;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.Set;
import javax.swing.*;
import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.table.DefaultTableCellRenderer;
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
	private HashSet<String> mutationFilterList= new HashSet<String>();
	
	private final int PANE_WIDTH= 380;
	private final int PANE_WIDTH_OFFSET= 20;
        final String CACHEFOLDER = System.getProperty("user.home")+"/.medsavant/plugins/cache/";
        final String OUTPUTDIR = "pathways_plugin/";
        final String GPMLFOLDER = "gpml_files/";
        final String HTMLFOLDER = "html_files/";
        private final int GENELIMIT = 2;
	private String[] mutationArray;
	private JPanel appView= new JPanel();
	private JPanel optionsPanel= new JPanel();
	private JScrollPane scrollPane= new JScrollPane();
	private JPanel resultsPanel= new JPanel();
        private JPanel pathwayInfoPanel = new JPanel();
	private IndividualSelector patientSelector= new IndividualSelector(true);
	private JButton choosePatientButton;
        private JButton performAnalysisButton;
        private JButton openPathwayButton;
        private JLabel pathwayListLabel;
        private JLabel tooFewGenesWarning;
	private String currentHospitalID;
	private String currentDNAID;
        private JComboBox pathwayList;
        private JComboBox multipleTestCorrectionList;
        private String[][] pathwayInfo;
        private JTabbedPane tabbedPane;
        private JPanel tablePanel;
        private JScrollPane ScrollPane;
        private JLabel pathwayImageLabel;
        private JLabel minPathwayGenesLabel;
        private JLabel maxPathwayGenesLabel;
        private JSlider minPathwayGenesSlider;
        private JSlider maxPathwayGenesSlider;
        private JSlider fdrSlider;
        private JLabel fdrLabel;
        private JLabel maxSliderLabel;
        private JLabel minSliderLabel;
        private JLabel fdrSliderLabel;
        private JFrame tableFrame;
        private String selectedPathwayTitle;
        private JPanel pathwayView;
        public static int BONFERRONI_INDEX = 0;
        public static int BENJAMINI_HOCHBERG_INDEX = 1;
        private final String[] multipleTestCorrections = {"Bonferroni","Benjamini-Hochberg (FDR)"};
        private int previouslySelectedTestCorrection;
        private SearchableTablePanel table;
        private Properties properties;
        private JLabel pathwayTitle;
        private JScrollPane pathwayDescriptionScrollPane;
        private JTextArea pathwayDescription;
        private JButton pngThumbnail;
        private JPanel infoPane;
        private JCheckBox lofCheckbox;
        private JCheckBox codingCheckbox;
        private JCheckBox nonsynonymousCheckbox;
        private static final String LOF_TEXT = "All loss of function mutations";
        private static final String CODING_TEXT = "All coding mutations";
        private static final String NONSYNONYMOUS_TEXT = "All non-synonymous mutations";
        private ArrayList<JCheckBox> mutationCheckBoxes;
        final private String PROPERTIESDIR = "properties/";
        final private String PROPERTIESFILENAME = "properties.xml";
        final private String DEFAULTPROPERTIESDIR = "medsavant/pathways/properties/default_properties.xml";
	private JButton maxSliderHelp;
        private JButton minSliderHelp;
        private JButton multipleTestCorrectionHelp;
        private JButton choosePatientHelp;
        private JButton exportGPML;
        private Font headerFont;
        private Font subHeaderFont;
        private Font boldFont;
        private final int headerFontSize = 24;
        private final int subHeaderFontSize = 20;
        private GregorianCalendar wikipathwaysDownloadVersion;
        private JScrollPane infoScrollPane;
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
		VariantType.ncRNA_INTRONIC.toString(), VariantType.ncRNA_SPLICING.toString()
	);
        private static final HashSet<String> LOF_MUTATIONS = new HashSet<String> (Arrays.asList(new String[] {
		VariantType.FS_DELETION.toString(), VariantType.FS_INSERTION.toString(),
		VariantType.FS_SUBSTITUTION.toString(), VariantType.FS_DUPLICATION.toString(),
		VariantType.SPLICING.toString(), VariantType.STOPGAIN.toString(),
                LOF_TEXT
        }));
        private static final HashSet<String> CODING_MUTATIONS= new HashSet<String> (Arrays.asList(new String[] {
		VariantType.MISSENSE.toString(), VariantType.SYNONYMOUS.toString(),
		VariantType.FS_DELETION.toString(), VariantType.FS_INSERTION.toString(),
		VariantType.FS_SUBSTITUTION.toString(), VariantType.FS_DUPLICATION.toString(),
		VariantType.NON_FS_DELETION .toString(), VariantType.NON_FS_INSERTION.toString(),
		VariantType.NON_FS_SUBSTITUTION.toString(), VariantType.NON_FS_DUPLICATION.toString(),
		VariantType.SPLICING.toString(), VariantType.STOPGAIN.toString(),
		VariantType.STOPLOSS.toString(), 
		VariantType.START_LOSS.toString(), CODING_TEXT
        }));
                
        private static final HashSet<String> NONSYNONYMOUS_MUTATIONS = new HashSet<String> (Arrays.asList(new String[] {
		VariantType.MISSENSE.toString(),
		VariantType.FS_DELETION.toString(), VariantType.FS_INSERTION.toString(),
		VariantType.FS_SUBSTITUTION.toString(), VariantType.FS_DUPLICATION.toString(),
		VariantType.NON_FS_DELETION .toString(), VariantType.NON_FS_INSERTION.toString(),
		VariantType.NON_FS_SUBSTITUTION.toString(), VariantType.NON_FS_DUPLICATION.toString(),
		VariantType.SPLICING.toString(), VariantType.STOPGAIN.toString(),
		VariantType.STOPLOSS.toString(), 
		VariantType.START_LOSS.toString(), NONSYNONYMOUS_TEXT
        }));
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
                //set properties - column widths & wikipathways download date
                this.getProperties();
                
                //set wikipathways download date
                wikipathwaysDownloadVersion = new GregorianCalendar(Integer.parseInt(properties.getProperty("wikipathwaysDownlaodYear")), Integer.parseInt(properties.getProperty("wikipathwaysDownlaodMonth")), Integer.parseInt(properties.getProperty("wikipathwaysDownlaodDay")));
                
                //initialize pathway info panel (contains pathway title, description, and png link-out)
                pathwayInfoPanel.setLayout(new MigLayout("fill"));
                pathwayTitle = new JLabel();
                
                Font font = pathwayTitle.getFont();
                boldFont = new Font(font.getFontName(), Font.BOLD, font.getSize());
                headerFont = new Font(font.getFontName(), Font.BOLD, headerFontSize);
                subHeaderFont = new Font(font.getFontName(), Font.BOLD, subHeaderFontSize);
                
                pathwayTitle.setFont(boldFont);
                pngThumbnail = new JButton();
                pngThumbnail.setToolTipText("Click to open pathway in browser");
                pngThumbnail.setHorizontalTextPosition(SwingConstants.CENTER);
                pngThumbnail.addActionListener(pngActionListener());
                pathwayDescriptionScrollPane = new JScrollPane();
                pathwayDescription = new JTextArea();
                pathwayDescription.setColumns(40);
                pathwayDescription.setLineWrap(true);
                pathwayDescription.setEditable(false);
                pathwayDescriptionScrollPane.setViewportView(pathwayDescription);
                exportGPML = new JButton("Export Pathway to .GPML for PathVisio");
                exportGPML.addActionListener(exportGPMLActionListener());
		// Create the options view
		optionsPanel.setLayout(new MigLayout("fillx"));
		optionsPanel.setMinimumSize(new Dimension(SIDE_PANE_WIDTH, 1));
		optionsPanel.setPreferredSize(new Dimension(SIDE_PANE_WIDTH, optionsPanel.getMaximumSize().height));
		optionsPanel.setBackground(LIGHTPURPLE);
		optionsPanel.setBorder(BorderFactory.createMatteBorder(0, 0, 0, 1, Color.LIGHT_GRAY));
		
		
		choosePatientButton= new JButton("Choose patient");
		choosePatientButton.addActionListener(choosePatientAction());
		choosePatientHelp = ViewUtil.getHelpButton("Choose Patient", 
				"Choose the variant call file that you would like the analysis to be run on. Variant call files can be uploaded using the VCF Upload app.");
		optionsPanel.add(choosePatientButton, "alignx center, split 2");
                
                optionsPanel.add(choosePatientHelp, "alignx center, wrap");
		
                
                mutationFilterGUI = this.mutationCheckboxPanel();
                mutationFilterGUI.setVisible(false);
                System.out.println("mutation width: "+mutationFilterGUI.getContentPaneWidth()+" height: "+mutationFilterGUI.getContentPaneHeight());
                System.out.println("Panel width: "+PANE_WIDTH+" offset: "+PANE_WIDTH_OFFSET);
                
                optionsPanel.add(mutationFilterGUI, "alignx center, wrap");
                optionsPanel.add(new JLabel(" "), "alignx center, wrap");
                
                previouslySelectedTestCorrection = 1;
		multipleTestCorrectionList = new JComboBox(multipleTestCorrections);
                multipleTestCorrectionList.setSelectedIndex(1);
                multipleTestCorrectionList.addActionListener(multipleTestCorrectionListener());
                multipleTestCorrectionList.setVisible(false);
                optionsPanel.add(multipleTestCorrectionList,"alignx center, split 2");
                multipleTestCorrectionHelp = ViewUtil.getHelpButton("Multiple Test Correction", 
				"The Bonferroni multiple test correction is more conservative, while the Benjamini Hochberg multiple test correction is less conservative. You may set a False Discovery Rate threshhold if you choose the Benjamini Hochberg.");
                multipleTestCorrectionHelp.setVisible(false);
                optionsPanel.add(multipleTestCorrectionHelp,"alignx center, wrap");
                
                fdrLabel = new JLabel("False discovery rate threshhold");
                fdrLabel.setVisible(false);
                optionsPanel.add(fdrLabel,"alignx center, wrap");
                fdrSliderLabel = new JLabel("0.05");
                fdrSliderLabel.setVisible(false);
                fdrSlider = new JSlider(0,50);
                fdrSlider.setValue(5);
                fdrSlider.addChangeListener(fdrSliderListener());
                fdrSlider.setVisible(false);
                optionsPanel.add(fdrSliderLabel,"alignx center, split 2");
                optionsPanel.add(fdrSlider,"alignx center, wrap");
                
                
                minSliderLabel = new JLabel("Minimum number of genes allowed in pathways");
                minSliderLabel.setVisible(false);
                optionsPanel.add(minSliderLabel,"alignx center, split 2");
                minSliderHelp = ViewUtil.getHelpButton("Minimum Pathway Size", 
				"Pathways with less genes than the threshhold you set will not be included in the pathway enrichment analysis.");
                minSliderHelp.setVisible(false);
                optionsPanel.add(minSliderHelp, "alignx center, wrap");
                
                minPathwayGenesLabel = new JLabel();
                minPathwayGenesSlider = new JSlider();
                minPathwayGenesSlider.addChangeListener(changeMinPathwayGenes());
                minPathwayGenesSlider.setVisible(false);
                optionsPanel.add(minPathwayGenesLabel,"alignx center, split 2");
                optionsPanel.add(minPathwayGenesSlider, "alignx center, wrap");
                
                maxSliderLabel = new JLabel("Maximum number of genes allowed in pathways");
                maxSliderLabel.setVisible(false);
                optionsPanel.add(maxSliderLabel,"alignx center, split 2");
                maxSliderHelp = ViewUtil.getHelpButton("Maximum Pathway Size", 
				"Pathways with more genes than the threshhold you set will not be included in the pathway enrichment analysis.");
                maxSliderHelp.setVisible(false);
                optionsPanel.add(maxSliderHelp,"alignx center, wrap");
                
                maxPathwayGenesLabel = new JLabel();
                maxPathwayGenesSlider = new JSlider();
                maxPathwayGenesSlider.addChangeListener(changeMaxPathwayGenes());
                maxPathwayGenesSlider.setVisible(false);
                optionsPanel.add(maxPathwayGenesLabel,"alignx center, split 2");
                optionsPanel.add(maxPathwayGenesSlider,"alignx center, wrap");
                
                optionsPanel.add(new JLabel(" "), "alignx center, wrap");
                
                performAnalysisButton = new JButton("Perform Hypergeometric Analysis");
                performAnalysisButton.addActionListener(analysisAction());
                performAnalysisButton.setVisible(false);
                optionsPanel.add(performAnalysisButton,"alignx center, wrap");
                
                tooFewGenesWarning = new JLabel();
                tooFewGenesWarning.setForeground(Color.RED);
                Font f = tooFewGenesWarning.getFont();
                tooFewGenesWarning.setFont(f.deriveFont(f.getStyle() | Font.BOLD));
                tooFewGenesWarning.setVisible(false);
                optionsPanel.add(tooFewGenesWarning,"alignx center, wrap");
                /*pathwayListLabel = new JLabel("");
                optionsPanel.add(pathwayListLabel, "alignx center, wrap");*/
                
                status= new JLabel();
		statusWheel= new ProgressWheel();
		statusWheel.setIndeterminate(true);
		// hide for now
		status.setVisible(false);
                
                /*pathwayList = new JComboBox(new String[0]);
                pathwayList.setPrototypeDisplayValue("this will be the width of the combo box");
                pathwayList.setVisible(false);
                optionsPanel.add(pathwayList,"alignx center, wrap");*/
                
                /*openPathwayButton = new JButton("Open selected pathway in browser");
                openPathwayButton.addActionListener(openPathwayAction());
                openPathwayButton.setVisible(false);
                optionsPanel.add(openPathwayButton, "alignx center, wrap");*/
                
		// Create the results view
		resultsPanel.setLayout(new MigLayout("fill"));
                resultsPanel.setBorder(BorderFactory.createMatteBorder(0, 0, 0, 0, Color.LIGHT_GRAY));
		resultsPanel.setBackground(Color.WHITE);
		
                resultPanelChoosePatientPrompt();
                
                infoPane = new JPanel();
                infoPane.setLayout(new MigLayout("insets 10 10 10 10, fillx"));
                infoScrollPane = new JScrollPane();
                infoScrollPane.setViewportView(infoPane);
                
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
	
        private ChangeListener changeMinPathwayGenes() {
            ChangeListener minChangeListener =  new ChangeListener() {
                @Override
		public void stateChanged(ChangeEvent e) {
                    //make analysis blank
                    resultPanelStartAnalysisPrompt();
                    //change min of max
                    maxPathwayGenesSlider.setMinimum(minPathwayGenesSlider.getValue());
                    //change label
                    minPathwayGenesLabel.setText(minPathwayGenesSlider.getValue()+"");
                }
            };
            return minChangeListener;
        }
        
        private ChangeListener changeMaxPathwayGenes() {
            ChangeListener maxChangeListener =  new ChangeListener() {
                @Override
		public void stateChanged(ChangeEvent e) {
                    //make analysis blank
                    resultPanelStartAnalysisPrompt();
                    //change max of min
                    minPathwayGenesSlider.setMaximum(maxPathwayGenesSlider.getValue());
                    //change label
                    maxPathwayGenesLabel.setText(maxPathwayGenesSlider.getValue()+"");
                }
                
            };
            return maxChangeListener;
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
                                        
                                        //process gene sets
                                        pathwayAnalysisObject.readWikiPathwayGeneSet();
                                        
                                        
					Set<String> selectedIndividuals= patientSelector.getHospitalIDsOfSelectedIndividuals();

					/* Once the user has made a patient hospital ID selection, get 
					 * the DNA ID so we can retrieve the patient's variants. */
					if (patientSelector.hasMadeSelection()) {
						currentHospitalID= patientSelector.getHospitalIDsOfSelectedIndividuals().iterator().next();
						String newDNAID= patientSelector.getDNAIDsOfSelectedIndividuals().iterator().next();

						if (newDNAID != null) {
							currentDNAID= newDNAID;
							choosePatientButton.setText(currentHospitalID);
                                                        mutationFilterGUI.setVisible(true);
                                                        multipleTestCorrectionList.setVisible(true);
                                                        
                                                        
                                                        resultPanelStartAnalysisPrompt();
                                                        
                                                        fdrSliderLabel.setVisible(true);
                                                        fdrLabel.setVisible(true);
                                                        fdrSlider.setVisible(true);
                                                        
                                                        int minPathwayGenes = pathwayAnalysisObject.getMinPathwayGenes();
                                                        int maxPathwayGenes = pathwayAnalysisObject.getMaxPathwayGenes();
                                                        System.out.println("MIN PATHWAY GENES: "+minPathwayGenes+" MAX PATHWAY GENES: "+maxPathwayGenes);
                                                        minPathwayGenesSlider.setMinimum(minPathwayGenes);
                                                        minPathwayGenesSlider.setMaximum(maxPathwayGenes);

                                                        if (minPathwayGenes < 10) {
                                                            minPathwayGenesSlider.setValue(10);
                                                        }
                                                        else {
                                                            minPathwayGenesSlider.setValue(minPathwayGenes);
                                                        }

                                                        minPathwayGenesLabel.setText(minPathwayGenesSlider.getValue()+"");

                                                        maxPathwayGenesSlider.setMinimum(minPathwayGenes);
                                                        maxPathwayGenesSlider.setMaximum(maxPathwayGenes);

                                                        if (maxPathwayGenes > 200) {
                                                            maxPathwayGenesSlider.setValue(200);
                                                        }
                                                        else {
                                                            maxPathwayGenesSlider.setValue(maxPathwayGenes);
                                                        }

                                                        maxPathwayGenesLabel.setText(maxPathwayGenesSlider.getValue()+"");

                                                        minPathwayGenesSlider.setVisible(true);
                                                        maxPathwayGenesSlider.setVisible(true);
                                                        minPathwayGenesLabel.setVisible(true);
                                                        maxPathwayGenesLabel.setVisible(true);
                                                        minSliderLabel.setVisible(true);
                                                        maxSliderLabel.setVisible(true);
                                                        optionsPanel.revalidate();
                                                        
                                                        //openPathwayButton.setVisible(false);
                                                        tooFewGenesWarning.setVisible(false);
                                                        optionsPanel.revalidate();
                                                        performAnalysisButton.setVisible(true);
                                                        optionsPanel.revalidate();
                                                        
                                                        multipleTestCorrectionHelp.setVisible(true);
                                                        minSliderHelp.setVisible(true);
                                                        maxSliderHelp.setVisible(true);
                                                        
						} else { // can't find this individual's DNA ID - may be a DB error
							errorDialog("Can't find a DNA ID for " + currentHospitalID);
						}
					}

                                        
				}
			};
		
		return outputAL;
	}
        private ChangeListener fdrSliderListener() {
            return new ChangeListener() {
                @Override
                public void stateChanged(ChangeEvent e) {
                    String labelText;
                    int sliderValue = fdrSlider.getValue();
                    if (sliderValue < 10) {
                        labelText = "0.0"+sliderValue;
                    }
                    else {
                        labelText = "0."+sliderValue;
                    }
                    fdrSliderLabel.setText(labelText);
                    optionsPanel.revalidate();
                }
            };
        }
        /**
	 * Update the variantPane with the set of variants.
	 */
	private ActionListener multipleTestCorrectionListener() {
            return new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    if (multipleTestCorrectionList.getSelectedIndex() == BENJAMINI_HOCHBERG_INDEX) {
                        fdrLabel.setVisible(true);
                        fdrSlider.setVisible(true);
                        fdrSliderLabel.setVisible(true);
                        optionsPanel.revalidate();
                    }
                    else {
                        fdrLabel.setVisible(false);
                        fdrSlider.setVisible(false);
                        fdrSliderLabel.setVisible(false);
                        optionsPanel.revalidate();
                    }
                    System.out.println(previouslySelectedTestCorrection + ", " + multipleTestCorrectionList.getSelectedIndex());
                    if (previouslySelectedTestCorrection!=multipleTestCorrectionList.getSelectedIndex()) {
                        resultPanelStartAnalysisPrompt();
                    }
                }
            };
        }
        private void updateResultsTable() {
		table.scrollSafeSelectAction(new Runnable() {
            @Override
            public void run() {
				
                if (table.getTable().getSelectedRow() != -1) {
                    SortableTable st= PathwaysPanel.this.table.getTable();
                    st.setAutoResizeMode(JTable.AUTO_RESIZE_SUBSEQUENT_COLUMNS);
                    int selectedIndex= st.getSelectedRow();
                    
                    String pathwayTitle = (String) st.getModel().getValueAt(selectedIndex, TestedPathway.PATHWAYNAMEINDEX);
                    
                    ImageIcon pngImage = pathwayAnalysisObject.getPathwayImage(pathwayTitle);
                    if (pngImage==null) {
                        pathwayImageLabel.setText("No image available for this pathway.");
                    }
                    else {
                        pathwayImageLabel.setIcon(pngImage);
                    }
                    selectedPathwayTitle = pathwayTitle;
                    pathwayInfoPanelLoadPathway(pathwayTitle, pngImage);
                    
                    
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
            //pathwayList.setVisible(false);
            status.setVisible(false);
            statusWheel.setVisible(false);
            //openPathwayButton.setVisible(false);
        }
        
        private void resultsPanelAddTabbedPane() {
            resultsPanel.removeAll();
            tabbedPane = ViewUtil.getMSTabedPane();
            tabbedPane.setBackground(Color.CYAN);
            tablePanel = new JPanel();
            tablePanel.setBackground(Color.PINK);
            tablePanel.setBorder(BorderFactory.createMatteBorder(0, 0, 0, 0, Color.LIGHT_GRAY));
            tablePanel.setLayout(new MigLayout("insets 0 0 0 0, fill"));
            tablePanel.add(table, "alignx center, gaptop 0, growx, growy, wrap 0");
            tablePanel.add(pathwayInfoPanel, "alignx center, growx, growy, gapbottom 0");
            tabbedPane.addTab("Results",tablePanel);
            pathwayImageLabel = new JLabel();
            tabbedPane.addTab("Information",infoScrollPane);
            tabbedPane.setToolTipTextAt(0,"P-values for the enrichment of each pathway");
            tabbedPane.setToolTipTextAt(1,"Information about analysis methods");
            resultsPanel.add(tabbedPane, "alignx center, growx, growy, wrap");
            pathwayInfoPanelPromptSelection();
            resultsPanel.revalidate();
            
            resultsPanel.repaint();
        }
        
        private void pathwayInfoPanelLoadPathway(String selectedPathwayTitle, ImageIcon pngImage) {
            pathwayInfoPanel.removeAll();
            pathwayTitle.setText(selectedPathwayTitle);
            pngThumbnail.setIcon(pngImage);
            pathwayDescription.setText(pathwayAnalysisObject.getDescription(selectedPathwayTitle));
            if (pathwayDescription.getText().equals("")) {
                pathwayInfoPanel.add(pathwayTitle,"alignx center, span 2");
                pathwayInfoPanel.add(pngThumbnail,"alignx center, growy, wrap");
            }
            else {
                pathwayInfoPanel.add(pathwayTitle,"alignx center, span 2");
                pathwayInfoPanel.add(pngThumbnail,"alignx center, growy, span 1 3, wrap");
                pathwayInfoPanel.add(pathwayDescriptionScrollPane,"alignx center, growx, growy, span 2, wrap");
                pathwayInfoPanel.add(exportGPML, "alignx center");
            }
            pathwayInfoPanel.revalidate();
        }
        
        private void pathwayInfoPanelPromptSelection() {
            pathwayInfoPanel.removeAll();
            pathwayInfoPanel.add(new JLabel("Please select pathway from table to see detailed pathway information."));
            pathwayInfoPanel.revalidate();
        }
        
        private ActionListener pngActionListener() {
            return new ActionListener () {
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
                                            pathwayAnalysisObject.generateGPMLFile(selectedPathwayTitle);
                                            pathwayAnalysisObject.generateHTMLFile(selectedPathwayTitle);
                                            //cancelLatch.countDown();
                                    } catch (Exception e) {
                                            errorDialog(e.getMessage());
                                            e.printStackTrace();
                                    }

                                    return null;
                            }

                            @Override
                            protected void showSuccess(Object t) {
                                //open link

                                String htmlPath = CACHEFOLDER+OUTPUTDIR+HTMLFOLDER;
                                String pathStart = "file://";
                                openLink(pathStart+htmlPath+pathwayAnalysisObject.getHTML(selectedPathwayTitle));
                                System.out.println("link open attempted");
                            }
                    };

                    // Execute thread
                    generateVisualizationThread.execute();  


                }
            };
        }
        
        private String fdrCutoffToString() {
            int fdr = fdrSlider.getValue();
            if (fdr < 10) {
                return "0.0"+fdr;
            }
            else {
                return "0."+fdr;
            }
        }
        
        private ActionListener exportGPMLActionListener() {
            return new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent ae) {
                    //open save dialog for user to save exported GPML file
                    JFileChooser save = new JFileChooser();
                    if (save.showSaveDialog(pathwayInfoPanel)==JFileChooser.APPROVE_OPTION) {
                        //make GPML
                        pathwayAnalysisObject.generateGPMLFile(selectedPathwayTitle);
                        
                        //copy to destination
                        try {
                         Path path = java.nio.file.Files.copy( 
                               new java.io.File(CACHEFOLDER+OUTPUTDIR+GPMLFOLDER+pathwayAnalysisObject.getGPML(selectedPathwayTitle)).toPath(), 
                               save.getSelectedFile().toPath(),
                               java.nio.file.StandardCopyOption.REPLACE_EXISTING,
                               java.nio.file.StandardCopyOption.COPY_ATTRIBUTES,
                               java.nio.file.LinkOption.NOFOLLOW_LINKS );
                        }
                        catch (Exception e) {
                            e.printStackTrace();
                        }
                    }
                }
            };
        }
        
        /**
	 * Action to perform when choose patient button is clicked.
	 * @return the ActionListener for this button
	 */
	/*private ActionListener openPathwayAction() {
		// create an anonymous class
		ActionListener outputAL= new ActionListener()
			{
                            @Override
                            public void actionPerformed(ActionEvent ae) {
                                //generate HTML files
                                
                                /* Background task. */
                                /*MedSavantWorker generateVisualizationThread= new MedSavantWorker<Object>(PathwaysPanel.class.getCanonicalName()) {			
                                        @Override
                                        protected Object doInBackground() {

                                                /* Create and perform a new analysis. Uses a CountDownLatch to 
                                                 * ensure that currentPGXAnalysis is initilized before I can
                                                 * do anything with it (for example, cancel it). */
                                                /*try {
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
	}*/
    
    /**
     * Open link in browser.
     * @param htmlFileLink String link file pathActionListener outputAL= new ActionListener()
			{
				@Override
				public void actionPerformed(ActionEvent ae) {
				//show status message and progress wheel
                                resultsPanelAnalyzing();
     */
    private static void openLink(String htmlFileLink) {
        try {
            //URI uri = new URI(htmlFileLink);
            
            //non functional attempt to use Cytoscape's catch-all open browser method
            
            
            //OpenBrowser.openURL(htmlFileLink);
            OpenBrowser openBrowser = new OpenBrowserImpl();
            openBrowser.openURL(htmlFileLink);
            //results in NoClassDefFoundError
            
            /*File file = new File(htmlFileLink);
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
            } else { System.out.println("desktop not supported - cannot open link in browsers"); }*/
            
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
	 * Action to perform when Perform Hypergeometric Analysis button is clicked.
	 * @return the ActionListener for this button
	 */
	private ActionListener analysisAction() {
		// create an anonymous class
		ActionListener outputAL= new ActionListener()
			{
				@Override
				public void actionPerformed(ActionEvent ae) {
				//show status message and progress wheel
                                resultsPanelAnalyzing();
                                    
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
                                        if (multipleTestCorrectionList.getSelectedIndex() == -1 || multipleTestCorrectionList.getSelectedIndex() == 0) {
                                            pathwayInfo = pathwayAnalysisObject.hypergeometricWithWikiPathways(currentDNAID, mutationFilterList, minPathwayGenesSlider.getValue(), maxPathwayGenesSlider.getValue(), multipleTestCorrectionList.getSelectedIndex());
                                        }
                                        else {
                                            double fdrCutoff = ((double) fdrSlider.getValue()) * 0.01;
                                            pathwayInfo = pathwayAnalysisObject.hypergeometricWithWikiPathways(currentDNAID, mutationFilterList, minPathwayGenesSlider.getValue(), maxPathwayGenesSlider.getValue(), multipleTestCorrectionList.getSelectedIndex(), fdrCutoff);
                                        }
					//cancelLatch.countDown();
				} catch (Exception e) {
					errorDialog(e.getMessage());
					e.printStackTrace();
				}

				return null;
			}

			@Override
			protected void showSuccess(Object t) {
                            table = pathwayAnalysisObject.getTableOutput();
                            resultsPanelAddTabbedPane();
                            choosePatientButton.setEnabled(true);
                            //make combo box with pathway titles for user to select from
                            //DELETE LATER - DO EVERYTHING FROM TABLE
                            //addPathwayTitlesToSelectionBox(pathwayInfo[0], pathwayInfo[3]);
                            //pathwayList.setVisible(true);
                            //add OK button
                            //openPathwayButton.setVisible(true);
                            
                            //create info tab
                            loadInfoPane();
                            
                            
                            int numGenes = pathwayAnalysisObject.genesInVariants();
                            if (numGenes < GENELIMIT) {
                                if (numGenes == 0) {
                                    tooFewGenesWarning.setText("No genes associated with variants.");
                                }
                                else if (numGenes == 1) {
                                    tooFewGenesWarning.setText("Warning: "+numGenes+" gene is too few for a quality analysis.");
                                }
                                else {
                                    tooFewGenesWarning.setText("Warning: "+numGenes+" genes is too few for a quality analysis.");
                                }
                                tooFewGenesWarning.setVisible(true);
                                optionsPanel.revalidate();
                            }
                            
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
        
        
        private void loadInfoPane() {
            JLabel largeTitleLabel = new JLabel("Pathway Analysis Notes and References");
            largeTitleLabel.setFont(headerFont);
            infoPane.add(largeTitleLabel, "alignx left, wrap");
            infoPane.add(new JLabel(" "), "alignx left, wrap");

            JLabel pathwayReferenceHeaderLabel = new JLabel("Source of Pathways");
            pathwayReferenceHeaderLabel.setFont(subHeaderFont);
            infoPane.add(pathwayReferenceHeaderLabel, "alignx left, wrap");
            
            SimpleDateFormat format = new SimpleDateFormat("MMM dd, yyyy");
            
            JLabel pathwayReferenceTextLabel = new JLabel("<html>Curated <i>Homo Sapiens</i> Pathway Collection from Wikipathways, downloaded on "+format.format(wikipathwaysDownloadVersion.getTime())+"</html>");
            infoPane.add(pathwayReferenceTextLabel, "alignx left, wrap");
            
            JButton wikipathwaysLinkButton = linkButton("Wikipathways Website", "http://www.wikipathways.org");
            infoPane.add(wikipathwaysLinkButton, "alignx left, wrap");
            infoPane.add(new JLabel(" "), "alignx left, wrap");
            //incl download dates.

            JLabel analysisTypeHeaderLabel = new JLabel("Analysis Type");
            analysisTypeHeaderLabel.setFont(subHeaderFont);
            infoPane.add(analysisTypeHeaderLabel, "alignx left, wrap");
            JTextArea analysisTypeText = new JTextArea("The analysis performed is the hypergeometric test. This determines the likelihood that, given a pathway size, a gene pool, and a list of variant genes, a pathway is enriched with genetic variants. Analysis is performed using the Hypergeometric class from the Java Statistical Classes package.");
            //figure out how to word that better
            analysisTypeText.setBackground(infoPane.getBackground());
            analysisTypeText.setLineWrap(true);
            analysisTypeText.setWrapStyleWord(true);
            infoPane.add(analysisTypeText, "alignx left, growx, wrap");
            JButton hypergeometricClassLink = linkButton("Hypergeometric Class from Java Statistical Classes","http://www.jsc.nildram.co.uk/api/jsc/distributions/Hypergeometric.html");
            infoPane.add(hypergeometricClassLink, "alignx left, wrap");
            JButton hypergeometricReferenceLink = linkButton("Wikipedia Hypergeometric Test Entry","http://en.wikipedia.org/wiki/Hypergeometric_distribution#Hypergeometric_test");
            infoPane.add(hypergeometricReferenceLink, "alignx left, wrap");
            
            
            infoPane.add(new JLabel(" "), "alignx left, wrap");

            JLabel analysisParametersHeaderLabel = new JLabel("Analysis Parameters");
            analysisParametersHeaderLabel.setFont(subHeaderFont);
            infoPane.add(analysisParametersHeaderLabel, "alignx left, wrap");
            infoPane.add(new JLabel("Only pathways containing between "+minPathwayGenesSlider.getValue() + " and "+ maxPathwayGenesSlider.getValue() + " genes were included in the analysis."), "alignx left, wrap");
            infoPane.add(new JLabel(" "), "alignx left, wrap");

            JLabel multipleTestCorrectionHeaderLabel = new JLabel("Multiple Test Correction");
            multipleTestCorrectionHeaderLabel.setFont(subHeaderFont);
            infoPane.add(multipleTestCorrectionHeaderLabel, "alignx left, wrap");
            if (multipleTestCorrectionList.getSelectedIndex() == BENJAMINI_HOCHBERG_INDEX) {
                infoPane.add(new JLabel("Benjamini-Hochberg multiple test correction, with a False Discovery Rate threshhold of "+fdrCutoffToString()), "alignx left, wrap");
            }
            else {
                infoPane.add(new JLabel("Bonferroni multiple test correction, corrected for "+pathwayAnalysisObject.getNumTestedPathways()+"total tests."), "alignx left, wrap");
            }
            infoPane.add(new JLabel(" "), "alignx left, wrap");
            //include FDR if benjamini hochberg

            JLabel visualizationHeaderLabel = new JLabel("Visualization");
            visualizationHeaderLabel.setFont(subHeaderFont);
            infoPane.add(visualizationHeaderLabel, "alignx left, wrap");
            JTextArea cytoscapeJSText = new JTextArea("Visualizations were created by converting the Wikipathways GPML format into pathways drawable by cytoscape.js.");
            cytoscapeJSText.setBackground(infoPane.getBackground());
            cytoscapeJSText.setLineWrap(true);
            cytoscapeJSText.setWrapStyleWord(true);
            infoPane.add(cytoscapeJSText, "alignx left, growx, wrap");
            JButton cytoscapeJSLink = linkButton("Cytoscape.js link","http://cytoscape.github.io/cytoscape.js/");
            infoPane.add(cytoscapeJSLink,"alignx left, wrap");
            infoPane.add(new JLabel(" "), "alignx left, wrap");
            //cytoscape.js

            JLabel analysisDetailsLabel = new JLabel("Analysis Details");
            analysisDetailsLabel.setFont(subHeaderFont);
            infoPane.add(analysisDetailsLabel,"alignx left, wrap");
            int genesUsed = pathwayAnalysisObject.numGenesInGeneSets();
            int genesNotUsed = pathwayAnalysisObject.numGenesNotInGeneSets();
            infoPane.add(new JLabel("Number of genes used in analysis: "+genesUsed),"alignx left, wrap");
            if (genesUsed > 0) {
                JScrollPane genesUsedScrollPane = new JScrollPane();
                JTextArea genesUsedText = new JTextArea(pathwayAnalysisObject.genesInGeneSetsText());
                genesUsedText.setEditable(false);
                genesUsedText.setRows(5);
                genesUsedScrollPane.setViewportView(genesUsedText);
                infoPane.add(genesUsedScrollPane,"alignx left, wrap");
            }
            infoPane.add(new JLabel("Number of genes not found in Wikipathways gene set: "+genesNotUsed),"alignx left, wrap");
            if (genesNotUsed > 0) {
                JScrollPane genesNotUsedScrollPane = new JScrollPane();
                JTextArea genesNotUsedText = new JTextArea(pathwayAnalysisObject.genesNotInGeneSetsText());
                genesNotUsedText.setEditable(false);
                genesNotUsedText.setRows(5);
                genesNotUsedScrollPane.setViewportView(genesNotUsedText);
                infoPane.add(genesNotUsedScrollPane,"alignx left, wrap");
            }
            infoPane.add(new JLabel(" "));
            infoPane.add(new JLabel("Number of pathways in Wikipathways gene set: "+pathwayAnalysisObject.getNumPathways()),"alignx left, wrap");
            JScrollPane allPathwaysScrollPane = new JScrollPane();
            JTextArea allPathwaysText = new JTextArea(pathwayAnalysisObject.allPathwaysText());
            allPathwaysText.setEditable(false);
            allPathwaysText.setRows(5);
            allPathwaysScrollPane.setViewportView(allPathwaysText);
            infoPane.add(allPathwaysScrollPane,"alignx left, wrap");
            infoPane.add(new JLabel("Number of pathways in Wikipathways gene set with associated genetic variants in this analysis: "+pathwayAnalysisObject.getNumTestedPathways()),"alignx left, wrap");
            if (pathwayAnalysisObject.getNumTestedPathways() > 0) {
                JScrollPane testedPathwaysScrollPane = new JScrollPane();
                JTextArea testedPathwaysText = new JTextArea(pathwayAnalysisObject.testedPathwaysText());
                testedPathwaysText.setEditable(false);
                testedPathwaysText.setRows(5);
                testedPathwaysScrollPane.setViewportView(testedPathwaysText);
                infoPane.add(testedPathwaysScrollPane,"alignx left, wrap");
            }
        }
        /**
	 * Fill combobox with pathway titles for selection
	 */
        /*private void addPathwayTitlesToSelectionBox(String[] pathwayTitles, String[] pValues) {
            for (int i = 0; i < pathwayTitles.length; i++) {
                if (pValues[i].equals("1")) {
                    pathwayList.addItem(pathwayTitles[i]);
                }
                else {
                    pathwayList.addItem(pValues[i]+" "+pathwayTitles[i]);
                }
                
            }
        }*/
        
	/**
	 * Create an error dialog and output the error to the log.
	 * @param errorMessage the error message to display.
	 */
	private void errorDialog(String errorMessage) {
		DialogUtils.displayError("Oops!", errorMessage);
		log.error("[" + this.getClass().getSimpleName() + "]: " + errorMessage);
	}
	
        private void resultsPanelAnalyzing () {
            resultsPanel.removeAll();
            resultsPanel.add(new JLabel("Analyzing pathway enrichment..."), "aligny bottom, alignx center, wrap");
            statusWheel.setVisible(true);
            resultsPanel.add(statusWheel, "aligny top, alignx center");
            resultsPanel.revalidate();
            resultsPanel.repaint();
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
                mutationCheckBoxes = new ArrayList<JCheckBox>();
                
                lofCheckbox = new JCheckBox(LOF_TEXT);
                lofCheckbox.addActionListener(lofActionListener());
                mutationCheckList.add(lofCheckbox, "wrap");
                mutationCheckBoxes.add(lofCheckbox);
                
                codingCheckbox = new JCheckBox(CODING_TEXT);
                codingCheckbox.addActionListener(codingActionListener());
                mutationCheckList.add(codingCheckbox, "wrap");
                mutationCheckBoxes.add(codingCheckbox);
                
                nonsynonymousCheckbox = new JCheckBox(NONSYNONYMOUS_TEXT);
                nonsynonymousCheckbox.addActionListener(nonsynonymousActionListener());
                mutationCheckList.add(nonsynonymousCheckbox, "wrap");
                mutationCheckBoxes.add(nonsynonymousCheckbox);
                
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
                                                if (mutationFilterList.equals(LOF_MUTATIONS)) {
                                                    lofCheckbox.setSelected(true);
                                                }
                                                else if (mutationFilterList.equals(CODING_MUTATIONS)) {
                                                    codingCheckbox.setSelected(true);
                                                }
                                                else if (mutationFilterList.equals(NONSYNONYMOUS_MUTATIONS)) {
                                                    nonsynonymousCheckbox.setSelected(true);
                                                }
                                                else {
                                                    lofCheckbox.setSelected(false);
                                                    codingCheckbox.setSelected(false);
                                                    nonsynonymousCheckbox.setSelected(false);
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
                        mutationCheckBoxes.add(currentCheckBox);
		}
		mutationCheckList.setBackground(LIGHTPURPLE);
		collapsibleMutation.setStyle(CollapsiblePane.PLAIN_STYLE);
		collapsibleMutation.setFocusPainted(false);
		collapsibleMutation.collapse(false);	
		
		collapsibleMutation.setMinimumSize(new Dimension(PANE_WIDTH - PANE_WIDTH_OFFSET, 0));
		collapsibleMutation.setContentPaneHeight(450);
                
                //add scroll
                JScrollPane scrollPane = new JScrollPane();
                scrollPane.setViewportView(mutationCheckList);
                collapsibleMutation.setContentPane(scrollPane);
		return collapsibleMutation;
	}
        
    private JButton linkButton(String text, String link) {
        final String hyperLink = link;
        JButton button = new JButton("<HTML><FONT color=\"#000099\"><U>"+text+"</U></FONT></HTML>");
            button.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    openLink(hyperLink);
                }
            });
        return button;
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
    private ActionListener lofActionListener() {
        return new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                        if (lofCheckbox.isSelected()) {
                            for (JCheckBox currentCheckBox : mutationCheckBoxes) {
                                if (LOF_MUTATIONS.contains(currentCheckBox.getText())) {
                                    currentCheckBox.setSelected(true);
                                    mutationFilterList.add(currentCheckBox.getText());
                                }
                                else {
                                    currentCheckBox.setSelected(false);
                                    mutationFilterList.remove(currentCheckBox.getText());
                                }
                            }
                        }
                }
        };
    }
    private ActionListener codingActionListener() {
        return new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                        if (codingCheckbox.isSelected()) {
                            for (JCheckBox currentCheckBox : mutationCheckBoxes) {
                                if (CODING_MUTATIONS.contains(currentCheckBox.getText())) {
                                    currentCheckBox.setSelected(true);
                                    mutationFilterList.add(currentCheckBox.getText());
                                }
                                else {
                                    currentCheckBox.setSelected(false);
                                    mutationFilterList.remove(currentCheckBox.getText());
                                }
                            }
                        }
                }
        };
    }
    private ActionListener nonsynonymousActionListener() {
        return new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                        if (nonsynonymousCheckbox.isSelected()) {
                            for (JCheckBox currentCheckBox : mutationCheckBoxes) {
                                if (NONSYNONYMOUS_MUTATIONS.contains(currentCheckBox.getText())) {
                                    currentCheckBox.setSelected(true);
                                    mutationFilterList.add(currentCheckBox.getText());
                                }
                                else {
                                    currentCheckBox.setSelected(false);
                                    mutationFilterList.remove(currentCheckBox.getText());
                                }
                            }
                        }
                }
        };
    }
    
    private void setTableColumnWidths() {
        SortableTable sortableTable = table.getTable();
        sortableTable.getColumnModel().getColumn(TestedPathway.PVALUEINDEX).setPreferredWidth(Integer.parseInt(properties.getProperty("pValueColumnWidth")));
        sortableTable.getColumnModel().getColumn(TestedPathway.PATHWAYNAMEINDEX).setPreferredWidth(Integer.parseInt(properties.getProperty("pathwayNameColumnWidth")));
        sortableTable.getColumnModel().getColumn(TestedPathway.GENESINDEX).setPreferredWidth(Integer.parseInt(properties.getProperty("genesColumnWidth")));
        sortableTable.getColumnModel().getColumn(TestedPathway.NUMGENESINDEX).setPreferredWidth(Integer.parseInt(properties.getProperty("numGenesColumnWidth")));
        
        sortableTable.getColumnModel().getColumn(TestedPathway.PVALUEINDEX).setCellRenderer(new DecimalFormatRenderer() );
    }
    
    private void saveNewProperties() {
        SortableTable sortableTable = table.getTable();
        properties.setProperty("pValueColumnWidth",sortableTable.getColumnModel().getColumn(TestedPathway.PVALUEINDEX).getWidth()+"");
        properties.setProperty("pathwayNameColumnWidth",sortableTable.getColumnModel().getColumn(TestedPathway.PATHWAYNAMEINDEX).getWidth()+"");
        properties.setProperty("genesColumnWidth",sortableTable.getColumnModel().getColumn(TestedPathway.GENESINDEX).getWidth()+"");
        properties.setProperty("numGenesColumnWidth",sortableTable.getColumnModel().getColumn(TestedPathway.NUMGENESINDEX).getWidth()+"");
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
	
        static class DecimalFormatRenderer extends DefaultTableCellRenderer {
		private static final DecimalFormat formatter = new DecimalFormat( "0.0000E0" );
 
		public Component getTableCellRendererComponent(
			JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
                        // First format the cell value as required

			value = formatter.format((Number)value);

                        // And pass it on to parent class 

			return super.getTableCellRendererComponent(
				table, value, isSelected, hasFocus, row, column );
		} 
	}
	
}




