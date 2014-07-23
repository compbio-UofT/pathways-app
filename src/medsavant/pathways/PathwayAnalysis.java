/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

//FOR TESTING PURPOSES, USE PATIENT 12-GM17206-Seqnom (many multiple genes)

package medsavant.pathways;

import com.healthmarketscience.sqlbuilder.BinaryCondition;
import com.healthmarketscience.sqlbuilder.ComboCondition;
import com.healthmarketscience.sqlbuilder.Condition;
import com.healthmarketscience.sqlbuilder.UnaryCondition;
import com.healthmarketscience.sqlbuilder.dbspec.basic.DbColumn;
import java.awt.Image;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.*;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.FileSystem;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.rmi.RemoteException;
import java.sql.SQLException;
import java.util.*;
import java.util.regex.*;
import javax.swing.ImageIcon;
import javax.xml.parsers.*;
import javax.xml.transform.*;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import javax.xml.xpath.*;
import jsc.distributions.Hypergeometric;
import org.springframework.core.io.Resource;
import org.springframework.core.io.support.PathMatchingResourcePatternResolver;
import org.ut.biolab.medsavant.MedSavantClient;
import org.ut.biolab.medsavant.client.project.ProjectController;
import org.ut.biolab.medsavant.client.reference.ReferenceController;
import org.ut.biolab.medsavant.client.util.DataRetriever;
import org.ut.biolab.medsavant.client.view.component.SearchableTablePanel;
import org.ut.biolab.medsavant.client.view.login.LoginController;
import org.ut.biolab.medsavant.shared.appdevapi.Variant;
import org.ut.biolab.medsavant.shared.appdevapi.VariantIterator;
import org.ut.biolab.medsavant.shared.db.TableSchema;
import org.ut.biolab.medsavant.shared.format.AnnotationFormat;
import org.ut.biolab.medsavant.shared.format.BasicVariantColumns;
import org.ut.biolab.medsavant.shared.model.SessionExpiredException;
import org.ut.biolab.medsavant.shared.serverapi.VariantManagerAdapter;
import org.w3c.dom.*;
import org.xml.sax.InputSource;
import org.apache.commons.io.FileUtils;
/**
 * This is a prototype for a MedSavant plugin that allows the user to
 *  perform a hypergeometric test on a list of variants and pathways
 *  and visualize pathways statistically likely to be affected by the variants
 *  in their browser.
 * @author ruthgrace
 */
public class PathwayAnalysis {
    //list of all pathways after hypergeometric test (pathways have p-values)
    private List<TestedPathway> testedPathways;
    //set of genes from variant list
    private Set<String> geneset;
    //all genes from genesets
    private Set<String> allgenes;
    //set of effects from variant list (is this ever used?)
    private Set<String> effectset;
    //list of all genesets from GeneSets or Wikipathways
    //  pathway titles are hashed to corresponding gene sets
    private HashMap genesets;
    //graphIDs of gene products in pathways hashed to their positions
    //  in the pathway visualization
    private HashMap<String,String> positionMap;
    //list of HTML links (local file links) for all pathways
    private String[] pathwayHtmlFileNames;
    //list of all pathway titles
    private String[] pathwayTitles;
    //list of all pathway file paths (GPML wikipathways format)
    private String[] pathwayGpmlFileNames;
    //Regular expressions for extracting variant information from VCF file
    final String GENE_REGEX="HGVS=([^:^(]*)";
    final String EFFECT_REGEX = "EFFECT=([^;]*)";
    //link to local test VCF files
    final String CACHEFOLDER = System.getProperty("user.home")+"/.medsavant/plugins/cache/";
    final String OUTPUTDIR = "pathways_plugin/";
    final String GPMLFOLDER = "gpml_files/";
    final String HTMLFOLDER = "html_files/";
    
    final String TESTFILE1 = "src/resources/annotatedVCF/WGS-001-03.gatk.snp.indel.jv.vcf";
    final String TESTFILE2 = "src/resources/annotatedVCF/WGS-002-03.gatk.snp.indel.jv.vcf";
    final String TESTFILE3 = "src/resources/annotatedVCF/WGS-003-03.gatk.snp.indel.jv.vcf";
    final String GENESETFOLDER = "src/genesets/";
    final String OUTPUTFILE = CACHEFOLDER+OUTPUTDIR+"enriched_pathways_and_Pvalues.txt";
    //final String PATHWAYOUTPUTFOLDER = "output/gpml_files/";
    final String WIKIPATHWAYSGMTFILE = "wikipathways.gmt";
    final String WIKIPATHWAYSFOLDER = "medsavant/pathways/wikipathwaysGPML/";
    final String PNGFOLDER = "medsavant/pathways/wikipathwaysPNG/";
    final String GENES_NOT_IN_GENESETS_FILE = CACHEFOLDER+OUTPUTDIR+"genes_not_in_genesets.txt";
    //final String PATHWAYSOUTPUTDIR = "output/gpml_files/";
    private static TableSchema ts= ProjectController.getInstance().getCurrentVariantTableSchema();
    private VariantManagerAdapter vma= MedSavantClient.VariantManager;
    private static final int DB_VARIANT_REQUEST_LIMIT= 500;
    private ArrayList<String> header;
    private List<PGXGene> pgxGenes= new LinkedList<PGXGene>();
    private HashMap<String,ImageIcon> pathwayLinks;
    private HashMap<String,String> pathwayGpmls;
    private int minPathwayGenes;
    private int maxPathwayGenes;
    private int multipleTestCorrection;
    private double fdrCutoff;
    private int minPathwayGenesFilter;
    private int maxPathwayGenesFilter;
    final public int NUM_TABLE_COLS=3;
    final double PVALUE_CUTOFF = 0.0000001; // you could turn this into something related to the number of pathways
    /**
     * constructor, initializes positionMap
     */
    public PathwayAnalysis() {
        header = new ArrayList<String>();
        header.add("P-values");
        header.add("Pathway Name");
        header.add("Genes associated\nwith patient variants\nin pathway");
        header.add("Genes\nin\npathway");
        positionMap = new HashMap<String,String>();
        //make output directories if they don't exist already
        File outputdir = new File(CACHEFOLDER+OUTPUTDIR+HTMLFOLDER);
        if (!outputdir.exists()) {
            outputdir.mkdir();
        }
        File gpmloutputdir = new File(CACHEFOLDER+OUTPUTDIR+GPMLFOLDER);
        if (!gpmloutputdir.exists()) {
            gpmloutputdir.mkdir();
        }
        File htmloutputdir = new File(CACHEFOLDER+OUTPUTDIR+HTMLFOLDER);
        if (!htmloutputdir.exists()) {
            htmloutputdir.mkdir();
        }
    }
    /**
     * Main method used to test a mode of analysis (Wikipathways or GeneSets)
     *  Wikipathways allows link outs to pathway visualizations.
     *  GeneSets only outputs a table of pathways and p-values
     * @param args 
     */
    /*public static void main(String[] args) {
        String WIKIPATHWAYSFOLDER = "src/resources/wikipathways/";
        String VCFFILE = "src/resources/annotatedVCF/WGS-001-03.gatk.snp.indel.jv.vcf";
        String GENESETFILE = "src/resources/genesets/wikipathways.gmt";
        PathwayAnalysis pathwayobject = new PathwayAnalysis();
        
        //run with wikipathways
        pathwayobject.hypergeometricWithWikiPathways("ERROR", new ArrayList<String>());
        
        //run with GeneSets
        //String GENESETFILE = "src/resources/genesets/Human_GO_AllPathways_no_GO_iea_symbol.gmt";
        //pathwayobject.hypergeometricWithGeneSets(GENESETFILE, VCFFILE);
        
    }*/
    /**
     * Run analysis with GeneSets (no visualizations), outputs table of
     *  pathway titles and p-values
     * @param GENESETFILE
     * @param VcfFile 
     */
    public void hypergeometricWithGeneSets(String GENESETFILE, String VcfFile) {
        //make output directories if they don't exist already
        File outputdir = new File(OUTPUTDIR);
        if (!outputdir.exists()) {
            outputdir.mkdir();
        }
        File gpmloutputdir = new File(CACHEFOLDER+OUTPUTDIR+GPMLFOLDER);
        if (!gpmloutputdir.exists()) {
            gpmloutputdir.mkdir();
        }
        File htmloutputdir = new File(CACHEFOLDER+OUTPUTDIR+HTMLFOLDER);
        if (!htmloutputdir.exists()) {
            htmloutputdir.mkdir();
        }
        //read in gene sets
        this.readGeneSet(GENESETFILE);
        //read in list of variants
        this.readVCF(VcfFile);
        //perform analysis
        this.hypergeometricTest();
        //output list of pathways and p-values
        this.outputEnrichedGeneList(OUTPUTFILE);
        //output basic stats
        System.out.println("Sample size: "+this.geneset.size()+", Population size: "+this.allgenes.size());
    }
    
    public String[][] hypergeometricWithWikiPathways(String currentDNAID, List<String> mutationTypes, int minPathwayGenesFilter, int maxPathwayGenesFilter) {
        return hypergeometricWithWikiPathways(currentDNAID, mutationTypes, minPathwayGenesFilter, maxPathwayGenesFilter, PathwaysPanel.BONFERRONI_INDEX);
    }
    /**
     * Run analysis with Wikipathways, allows link-out visualizations of 
     *  pathways suspected to be affected
     * @param VcfFile
     * @param pathwayfolder 
     */
    public String[][] hypergeometricWithWikiPathways(String currentDNAID, List<String> mutationTypes, int minPathwayGenesFilter, int maxPathwayGenesFilter, int multipleTestCorrection) {
        this.multipleTestCorrection = multipleTestCorrection;
        this.minPathwayGenesFilter = minPathwayGenesFilter;
        this.maxPathwayGenesFilter = maxPathwayGenesFilter;
        
        try {
            //this.wikipathways2GMT();
            this.queryVariants(currentDNAID, mutationTypes);
            //perform analysis
            this.hypergeometricWikiPathwaysTest();
            //output pathways and p-values
            this.outputEnrichedGeneList(OUTPUTFILE);
            //make all non-inhibitor non-metabolite non-mutated gene products black
            //this.makeGPMLnodesblack(PATHWAYOUTPUTFOLDER);
            /*String[][] pathwayInfo = new String[4][pathwayTitles.length];
            pathwayInfo[0] = pathwayTitles;
            pathwayInfo[1] = pathwayHtmlFileNames;
            pathwayInfo[2] = pathwayGpmlFileNames;
            pathwayInfo[3] = pathwayPValues;*/
            this.assignPNGFiles();
            System.out.println("Sample size: "+this.geneset.size()+", Population size: "+this.allgenes.size());
            
            return TestedPathway.getPathwayInfoArray(testedPathways);
        }
        catch (Exception e) {
            e.printStackTrace();
            return null;
        }
    }
    
    /**
     * Run analysis with Wikipathways, allows link-out visualizations of 
     *  pathways suspected to be affected
     * @param VcfFile
     * @param pathwayfolder 
     */
    public String[][] hypergeometricWithWikiPathways(String currentDNAID, List<String> mutationTypes, int minPathwayGenesFilter, int maxPathwayGenesFilter, int multipleTestCorrection, double fdr) {
        this.multipleTestCorrection = multipleTestCorrection;
        this.minPathwayGenesFilter = minPathwayGenesFilter;
        this.maxPathwayGenesFilter = maxPathwayGenesFilter;
        this.fdrCutoff = fdr;
        try {
            //this.wikipathways2GMT();
            this.queryVariants(currentDNAID, mutationTypes);
            //perform analysis
            this.hypergeometricWikiPathwaysTest();
            //output pathways and p-values
            this.outputEnrichedGeneList(OUTPUTFILE);
            //make all non-inhibitor non-metabolite non-mutated gene products black
            //this.makeGPMLnodesblack(PATHWAYOUTPUTFOLDER);
            /*String[][] pathwayInfo = new String[4][pathwayTitles.length];
            pathwayInfo[0] = pathwayTitles;
            pathwayInfo[1] = pathwayHtmlFileNames;
            pathwayInfo[2] = pathwayGpmlFileNames;
            pathwayInfo[3] = pathwayPValues;*/
            this.assignPNGFiles();
            System.out.println("Sample size: "+this.geneset.size()+", Population size: "+this.allgenes.size());
            
            return TestedPathway.getPathwayInfoArray(testedPathways);
        }
        catch (Exception e) {
            e.printStackTrace();
            return null;
        }
    }
    
    
    private void assignPNGFiles() {
        this.pathwayLinks = new HashMap<String,ImageIcon>();
        String path;
        for (int i = 0; i < pathwayTitles.length; i++) {
            try {
                path = pathwayHtmlFileNames[i].replaceFirst("html","png");
                pathwayLinks.put(pathwayTitles[i],new ImageIcon(getClass().getClassLoader().getResource(PNGFOLDER+path)));
            }
            catch (NullPointerException e) {
                System.out.println("no resource at path: "+PNGFOLDER+pathwayHtmlFileNames[i].replaceFirst("html","png"));
                
            }
        }
    }
    
    /**
     * process gene sets (derived from Wikipathways GPML files) so that
     *  they are stored in this object
     * @param filename 
     */
    public void readWikiPathwayGeneSet() {
        //store genesets as hashmap of hashset, hashmap is the pathway name.
        //store genes as hashmap of ArrayList filled with pathways.
        genesets = new HashMap<String,HashMap<String,Integer>>();
        minPathwayGenes = Integer.MAX_VALUE;
        maxPathwayGenes = Integer.MIN_VALUE;
        allgenes = new HashSet<String>();
        pathwayGpmls = new HashMap<String,String>();
        String[] linecomponents = new String[0];
        try {
            //BufferedReader br = new BufferedRßeader(new FileReader(this.GENESETFOLDER+filename));
            BufferedReader br = new BufferedReader( new InputStreamReader(this.getClass().getResourceAsStream("/medsavant/pathways/geneSet/wikipathways.gmt")));
            String line = br.readLine();
            String pathwayName;
            HashMap<String,Integer> pathwaygenes;
            ArrayList<String> names = new ArrayList<String>();
            ArrayList<String> gpmls = new ArrayList<String>();
            ArrayList<String> htmls = new ArrayList<String>();
            while (line != null) {
                line = line.trim();
                //tab separated
                //caps name % humancyc % abbreviation, human readable name, gene symbols
                linecomponents = line.split("\t");
                if (linecomponents.length >2) {
                    
                    pathwayName = linecomponents[1].split("\\|")[0];
                    names.add(pathwayName);
                    gpmls.add(linecomponents[0]);
                    htmls.add(linecomponents[0].replaceFirst("gpml","html"));
                    pathwaygenes = new HashMap<String,Integer>();
                    int counter = 0;
                    for (int i = 2; i < linecomponents.length; i++) {
                        pathwaygenes.put(linecomponents[i],counter++);
                        allgenes.add(linecomponents[i]);
                    }
                    this.genesets.put(pathwayName,pathwaygenes);
                    pathwayGpmls.put(pathwayName, linecomponents[0]);
                    if (pathwaygenes.size() < minPathwayGenes && pathwaygenes.size()!=0){
                        minPathwayGenes = pathwaygenes.size();
                    }
                    if (pathwaygenes.size() > maxPathwayGenes){
                        maxPathwayGenes = pathwaygenes.size();
                    }
                }
                else {
                    System.out.println("unused genesets line: "+line);
                }
                line = br.readLine();
            }   
            pathwayTitles = new String[1];
            this.pathwayTitles = names.toArray(pathwayTitles);
            this.pathwayGpmlFileNames = new String[1];
            this.pathwayGpmlFileNames = gpmls.toArray(pathwayGpmlFileNames);
            this.pathwayHtmlFileNames = new String[1];
            this.pathwayHtmlFileNames = htmls.toArray(pathwayHtmlFileNames);
            System.out.println("calculated min and max pathway genes are: MIN - " +minPathwayGenes+" MAX - "+maxPathwayGenes);
        }
        catch (Exception e) {
            for (int counter = 0; counter < linecomponents.length; counter++) {
                System.out.println(linecomponents[counter]);
            }
            
            e.printStackTrace();
        }
    }
    
    /** 
	 * Searchable table output for development testing. 
	 * @param selectedViewColumns Columns preselected for SearchableTablePanel output
	 */
	public SearchableTablePanel getTableOutput() {
            int[] selectedViewColumns = {0,1,2,3};
		DataRetriever<Object[]> dr= new DataRetriever<Object[]>() {
			@Override
			public List<Object[]> retrieve(int start, int limit) throws Exception {            
				//return allVariants;
                            if (multipleTestCorrection == PathwaysPanel.BONFERRONI_INDEX) {
                                return TestedPathway.convertToObjectListBonferroni(testedPathways);
                            }
                            else if (multipleTestCorrection == PathwaysPanel.BENJAMINI_HOCHBERG_INDEX) {
                                return TestedPathway.convertToObjectListBH(testedPathways, fdrCutoff);
                            }
                            else {
                                System.out.println("INVALID MULTIPLE TEST CORRECTION");
                            }
                            return null;
			}

			@Override
			public int getTotalNum() {
				//return allVariants.size();
                            return testedPathways.size();
			}

			@Override
			public void retrievalComplete() {
			}
		};
		
                if (multipleTestCorrection == 1) {
                    header.set(0, "Q-values");
                }
                
		Class[] STRING_ONLY_COLUMN_CLASSES= new Class[header.size()];
		for (int i= 0; i != STRING_ONLY_COLUMN_CLASSES.length; ++i)
			STRING_ONLY_COLUMN_CLASSES[i]= String.class; // FOR NOW ONLY CALLING THESE STRINGS
		
		SearchableTablePanel t;
		
		// if the selected columns use incorrect/outdated indices, default to all columns
		try {
			
			if (selectedViewColumns == null) {
				t= new SearchableTablePanel("Results", header.toArray(new String[header.size()]), 
					STRING_ONLY_COLUMN_CLASSES, new int[0], true, true, Integer.MAX_VALUE,
					false, SearchableTablePanel.TableSelectionType.ROW, Integer.MAX_VALUE, dr);
			} else {
				t= new SearchableTablePanel("Results", header.toArray(new String[header.size()]), 
					STRING_ONLY_COLUMN_CLASSES, new int[0], 
					true, true, Integer.MAX_VALUE, false, 
					SearchableTablePanel.TableSelectionType.ROW, Integer.MAX_VALUE, dr);
			}
		} catch (Exception e) {
			t= new SearchableTablePanel("Results", header.toArray(new String[header.size()]), 
				STRING_ONLY_COLUMN_CLASSES, new int[0], true, true, Integer.MAX_VALUE,
				false, SearchableTablePanel.TableSelectionType.ROW, Integer.MAX_VALUE, dr);
		}
		
                
		t.setResizeOff();
		t.setExportButtonVisible(true);
		t.setExportButtonEnabled(true);
		t.setHelpButtonVisible(false);
		//t.setChooseColumnsButtonVisible(false);
		t.forceRefreshData(); // without this, the table is empty with just a header
		
		return t;
	}
    
    public int getMaxPathwayGenes() {
        return this.maxPathwayGenes;
    }
    
    public int getMinPathwayGenes() {
        return this.minPathwayGenes;
    }
    
    /**
     * process gene sets (derived from GeneSets gmt files) so that
     *  they are stored in this object
     * @param filename 
     */
    public void readGeneSet(String filename) {
        //store genesets as hashmap of hashset, hashmap is the pathway name.
        //store genes as hashmap of ArrayList filled with pathways.
        genesets = new HashMap<String,ArrayList>();
        allgenes = new HashSet<String>();
        String[] linecomponents;
        try {
            BufferedReader br = new BufferedReader(new FileReader(filename));
            String line = br.readLine();
            String pathwayname;
            ArrayList pathwaygenes;
            while (line != null) {
                line = line.trim();
                //tab separated
                //caps name % humancyc % abbreviation, human readable name, gene symbols
                linecomponents = line.split("\t");
                pathwayname = linecomponents[0] + "\t" + linecomponents[1].split("\\|")[0];
                pathwaygenes = new ArrayList<String>();
                for (int i = 2; i < linecomponents.length; i++) {
                    pathwaygenes.add(linecomponents[i]);
                    allgenes.add(linecomponents[i]);
                }
                this.genesets.put(pathwayname,pathwaygenes);
                line = br.readLine();
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    /**
     * make all nodes black in GPML files
     * @param folderpath 
     */
/*    public void makeGPMLnodesblack(String gpmlFileName) {
        try {
            final File folder = new File(folderpath);
            for (final File fileEntry : folder.listFiles()) {
                Path path = Paths.get(folderpath+fileEntry.getName());
                Charset charset = StandardCharsets.UTF_8;
                String content = new String(Files.readAllBytes(path), charset);
                //REPLACE COLOR REGEX WITH COLOR BLACK
                content = content.replaceAll("foo", "bar");
                Files.write(path, content.getBytes(charset));
                
                fileEntry.getName();
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }*/
    /**
     * Generate HTML visualization of pathway, rendered by converting
     *  Wikipathways GPML into javascript displayable by cytoscape.js
     *  with mutated genes colored green with thick node outlines
     *  metabolites colored blue, and inhibitors colored red.
     */
    public void generateHTMLFile(String gpmlPath) {
        String outputFolder = CACHEFOLDER+OUTPUTDIR;
        
        //extract nodes (Label = gene symbol, ID = something unique, use position, color = black unless in tested pathways)
        try {
        
            Document doc;
            XPath xpath;
            PrintWriter writer;
            String xpathexpression, geneSymbol;
            NodeList nodes;
            Node node, graphics,point;
            String[] javascriptNodes=new String[1];
            String[] javascriptEdges=new String[1];
            
            NamedNodeMap attributes;
            String htmlFilePath, pathwayName, pathwayDescription;
            
            
            int counter = 0;
            HashMap<String,String> nonGeneProductNodes;
            HashMap<String,String> geneProductNodes;
            
                doc = DocumentBuilderFactory.newInstance().newDocumentBuilder().parse(new FileInputStream(outputFolder+GPMLFOLDER+gpmlPath));
                xpath = XPathFactory.newInstance().newXPath();
                nonGeneProductNodes = new HashMap<String,String>();
                geneProductNodes = new HashMap<String,String>();
                xpathexpression = "/Pathway";
                nodes = (NodeList) xpath.evaluate(xpathexpression,doc,XPathConstants.NODESET);
                
                if (nodes!=null && nodes.getLength() >0) {
                    attributes = nodes.item(0).getAttributes();
                    pathwayName=attributes.getNamedItem("Name").getTextContent();
                }
                else {
                    System.out.println("could not find pathway name of "+outputFolder+GPMLFOLDER+gpmlPath);
                    pathwayName="No Pathway Name";
                }
                htmlFilePath = outputFolder+GPMLFOLDER+gpmlPath.replaceFirst("\\.gpml","\\.html");
                
                xpathexpression = "/Pathway/Comment";
                nodes = (NodeList) xpath.evaluate(xpathexpression,doc,XPathConstants.NODESET);
                if (nodes.getLength() >0 && nodes.item(0).getNodeType()==Node.ELEMENT_NODE){
                    int commentNodeCounter = 0;
                    Node commentNode = nodes.item(commentNodeCounter);
                    while (commentNode!=null && commentNode.getAttributes().getNamedItem("Source")!= null && !commentNode.getAttributes().getNamedItem("Source").getNodeValue().equals("WikiPathways-description")){
                        commentNode = nodes.item(++commentNodeCounter);
                    }
                    if (commentNode!=null) {
                        pathwayDescription = commentNode.getTextContent();
                    }
                    else {
                        pathwayDescription = "";
                    }
                }
                else {
                    pathwayDescription = "";
                }
                
                counter++;
                
                xpathexpression = "/Pathway/DataNode";
                nodes = (NodeList) xpath.evaluate(xpathexpression,doc,XPathConstants.NODESET);
                //gene products
                if (nodes.getLength() >0 && nodes.item(0).getNodeType()==Node.ELEMENT_NODE){
                    javascriptNodes = new String[nodes.getLength()];
                    for (int i = 0; i < javascriptNodes.length; i++) {
                        node = nodes.item(i);
                        javascriptNodes[i] = this.processNode(node, nonGeneProductNodes, geneProductNodes, this.positionMap);
                    }
                }
                
                //metabolites, etc (everything but gene products)
                xpathexpression = "/Pathway/Label";
                nodes = (NodeList) xpath.evaluate(xpathexpression,doc,XPathConstants.NODESET);
                 if (nodes.getLength() >0 && nodes.item(0).getNodeType()==Node.ELEMENT_NODE){
                     for (int i = 0; i < javascriptNodes.length; i++) {
                        node = nodes.item(i);
                        if (node == null || !node.hasChildNodes()) {
                            continue;
                        }
                        this.processNode(node, nonGeneProductNodes, geneProductNodes, this.positionMap);
                     }
                     //remove trailing comma
                 }
                
                //interactions
                xpathexpression = "/Pathway/Interaction";
                nodes = (NodeList) xpath.evaluate(xpathexpression,doc,XPathConstants.NODESET);
                //get nodes
                if (nodes.getLength() >0 && nodes.item(0).getNodeType()==Node.ELEMENT_NODE){
                    javascriptEdges = new String[nodes.getLength()];
                    for (int i = 0; i < javascriptEdges.length; i++) {
                        node = nodes.item(i);
                        graphics = node.getFirstChild();
                        while (graphics!=null && !graphics.getNodeName().equals("Graphics")) {
                            graphics = graphics.getNextSibling();
                        }
                        if (graphics == null ) {
                            System.out.println("edge did not have graphics component");
                            continue;
                        }
                        javascriptEdges[i] = this.processEdge(graphics, nonGeneProductNodes, geneProductNodes);
                    }
                }
                
                //shapes
                xpathexpression = "/Pathway/Shape";
                nodes = (NodeList) xpath.evaluate(xpathexpression,doc,XPathConstants.NODESET);
                 if (nodes.getLength() >0 && nodes.item(0).getNodeType()==Node.ELEMENT_NODE){
                     for (int i = 0; i < javascriptNodes.length; i++) {
                        node = nodes.item(i);
                        if (node == null || !node.hasChildNodes()) {
                            continue;
                        }
                        this.processNode(node, nonGeneProductNodes, geneProductNodes, this.positionMap);
                     }
                 }
                
                writeJS(pathwayName, pathwayDescription, outputFolder+HTMLFOLDER, gpmlPath.replaceFirst("\\.gpml","\\.js"), javascriptNodes, javascriptEdges, nonGeneProductNodes);
                writeHTML(gpmlPath, outputFolder+HTMLFOLDER, pathwayName, pathwayDescription);
            
            this.writeCSS(outputFolder+HTMLFOLDER+"cytoscape_javascript_prototype.css");
            this.writeShowMore(outputFolder+HTMLFOLDER);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
    }
    /**
     * output javascript file for HTML link out
     * @param pathwayName
     * @param pathwayDescription
     * @param folder
     * @param fileName
     * @param javascriptNodes
     * @param javascriptEdges
     * @param nonGeneProductNodes 
     */
    private void writeJS(String pathwayName, String pathwayDescription, String folder, String fileName, String[] javascriptNodes, String[] javascriptEdges, HashMap<String,String> nonGeneProductNodes) {
        try {
            String jsfilebeginning = "$('#cy').cytoscape({\n  layout: {\n  name: 'preset',\n  positions: {";
            String jsfilestyle = "}\n},\nstyle: cytoscape.stylesheet()\n    .selector('node')\n      .css({\n        'content': 'data(name)',\n        'shape': 'rectangle',\n        'text-valign': 'center',\n        'background-color': 'white',\n        'background-opacity': 1,\n        'color': 'black',\n        'font-size': 10,\n        'text-outline-width': 0,\n        'border-width': 2,\n        'border-color': 'black',\n        'border-opacity': 1,\n        'text-outline-color': '#888',\n        'z-index': 2\n      })\n      .selector('$node > node')\n      .css({\n        'padding-top': '4px',\n        'padding-left': '4px',\n        'padding-bottom': '4px',\n        'padding-right': '4px',\n        'z-index': 0\n      })\n      .selector('node[isInhibitor]')\n      .css({\n        'color': 'red',\n        'border-color': 'red'\n      })\n      .selector('node[isMetabolite]')\n      .css({\n        'color': 'blue',\n        'border-color': 'blue'\n      })\n      .selector('node[placementNode]')\n      .css({\n        'width': 1,\n        'height': 1,\n        'background-opacity': 1,\n        'background-color': 'black',\n        'border-opacity': 0\n      })\n      .selector('node[width]')\n      .css({\n        'width': 'data(width)',\n        'height': 'data(height)',\n      })\n      .selector('node[shape]')\n      .css({\n        'background-opacity': 0,\n        'background-color': 'white',\n        'shape': 'data(shape)',\n        'z-index': 0\n      })\n      .selector('node > $node')\n      .css({\n        'z-index': 0\n      })\n      .selector('node[isLabel]')\n      .css({\n      	'background-opacity': 0,\n      	'border-opacity': 0,\n        'z-index': 1\n      })\n      .selector('node[color]')\n      .css({\n        'border-color': 'data(color)'\n      })\n      .selector('node[isMutated]')\n      .css({\n        'border-width': 4,\n        'border-color': 'green'\n      })\n      .selector('edge')\n      .css({\n        'target-arrow-shape': 'triangle',\n        'z-index': 3,\n        'width': 2,\n        'line-color': 'black',\n        'target-arrow-color': 'black'\n      })\n    .selector('edge[noArrowHead]')\n      .css({\n        'target-arrow-shape': 'none'\n      })\n    .selector('edge[teeArrowHead]')\n      .css({\n        'target-arrow-shape': 'tee',\n        'line-color': 'red',\n        'target-arrow-color': 'red'\n      })\n    .selector('edge[veeArrowHead]')\n      .css({\n        'target-arrow-shape': 'vee'\n      })\n    .selector('edge[color]')\n      .css({\n        'line-color': 'data(color)'\n      })\n    .selector('edge[circularArrowHead]')\n      .css({\n        'target-arrow-shape': 'circle',\n        'target-arrow-fill': 'hollow'\n      })\n    .selector(':selected')\n      .css({\n        'background-color': 'black',\n        'line-color': 'black',\n        'target-arrow-color': 'black',\n        'source-arrow-color': 'black'\n      })\n    .selector('.faded')\n      .css({\n        'background-opacity': 0.25,\n        'text-opacity': 0\n      })\n    .selector('node#supergrandparent')\n      .css({\n        'background-opacity': 0,\n        'border-opacity': 0,\n        'border-color': 'white',\n        'z-index': -1\n      })\n    .selector('node#superparent')\n      .css({\n        'background-opacity': 0,\n        'border-opacity': 0,\n        'border-color': 'white',\n        'z-index': -1\n      }),\nelements: {\n    nodes: [\n    { data: { id: 'supergrandparent' } },\n    { data: { id: 'superparent', parent: 'supergrandparent'} },";
            String jsfilemiddle = "    ],\n    edges: [";
            String jsfileend = "    ]\n  },\n  \n  ready: function(){\n    window.cy = this;\n    \n    // giddy up...\n    \n    cy.elements().unselectify();\n    \n    cy.on('tap', 'node', function(e){\n      var node = e.cyTarget; \n      var neighborhood = node.neighborhood().add(node);\n      \n      cy.elements().addClass('faded');\n      neighborhood.removeClass('faded');\n    });\n    \n    cy.on('tap', function(e){\n      if( e.cyTarget === cy ){\n        cy.elements().removeClass('faded');\n      }\n    });\n  }\n});";
            PrintWriter writer = new PrintWriter(folder+fileName, "UTF-8");
            writer.println(jsfilebeginning);
            //print node position map
            Iterator positions = positionMap.values().iterator();
            if (positions.hasNext()) {
                writer.println(positions.next());
            }
            while (positions.hasNext()) {
                writer.println(", "+positions.next());
            }
            writer.println(jsfilestyle);
            writer.print(javascriptNodes[0]);
            for (int i = 1; i < javascriptNodes.length; i++) {
                if (!javascriptNodes[i].equals("")) {
                    writer.print(",\n"+javascriptNodes[i]);
                }
            }
            
            Iterator<String> it = nonGeneProductNodes.values().iterator();
            while (it.hasNext()) {
                writer.print(",\n"+it.next());
            }
            writer.println();
            writer.println(jsfilemiddle);
            
            //remove trailing comma
            if (javascriptEdges[javascriptEdges.length-1].substring(javascriptEdges[javascriptEdges.length-1].length() - 1).equals(",")) {
                javascriptEdges[javascriptEdges.length-1] = javascriptEdges[javascriptEdges.length-1].substring(0,javascriptEdges[javascriptEdges.length-1].length() - 1);
                
            }
            
            for (int i = 0; i < javascriptEdges.length; i++) {
                writer.println(javascriptEdges[i]);
            }
            writer.println(jsfileend);
            writer.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
    }
    
    public int getNumGenesInPathway(String pathwayName) {
        return ( (HashMap)genesets.get(pathwayName) ).size();
    }
    private void writeShowMore(String folder) {
        try {
            PrintWriter writer = new PrintWriter(folder+"showmore.js");
            writer.print("$(document).ready(function() {\n    $(\".show-more a\").each(function() {\n    var $link = $(this);\n    var $content = $link.parent().prev(\"div.text-content\");\n\n    console.log($link);\n\n    var visibleHeight = $content[0].clientHeight;\n    var actualHide = $content[0].scrollHeight - 1;\n\n    console.log(actualHide);\n    console.log(visibleHeight);\n\n    if (actualHide > visibleHeight) {\n        $link.show();\n    } else {\n        $link.hide();\n    }\n});\n\n$(\".show-more a\").on(\"click\", function() {\n    var $link = $(this);\n    var $content = $link.parent().prev(\"div.text-content\");\n    var linkText = $link.text();\nconsole.log(\"click!\");\n    $content.toggleClass(\"short-text, full-text\");\n\n    $link.text(getShowLinkText(linkText));\n\n    return false;\n});\n\nfunction getShowLinkText(currentText) {\n    var newText = '';\n\n    if (currentText.toUpperCase() === \"SHOW MORE\") {\n        newText = \"Show less\";\n    } else {\n        newText = \"Show more\";\n    }\n\n    return newText;\n}\n});\n\n\n");
            writer.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    /**
     * output HTML file for pathway visualization in browser
     * @param gpmlFileName
     * @param folder
     * @param pathwayName
     * @param pathwayDescription 
     */
    private void writeHTML(String gpmlFileName, String folder, String pathwayName, String pathwayDescription) {
        try {
            PrintWriter writer = new PrintWriter(folder+gpmlFileName.replaceFirst("\\.gpml","\\.html"));
            
            
            writer.print("<!DOCTYPE html>\n<html>\n<head>\n<meta name=\"description\" content=\"[Pathway Display]\" />\n<script src=\"http://ajax.googleapis.com/ajax/libs/jquery/1/jquery.min.js\"></script>\n<meta charset=utf-8 />\n<title>");
            writer.print(pathwayName);
            writer.print("</title>\n  <script src=\"http://cytoscape.github.io/cytoscape.js/api/cytoscape.js-latest/cytoscape.min.js\"></script>\n<script src=\"showmore.js\"></script>\n  <link rel=\"stylesheet\" href=\"./cytoscape_javascript_prototype.css\" />\n\n<div class=\"text-container\">\n<h1>");
            writer.print(pathwayName+"</h1><div class=\"text-content short-text\">");
            writer.println(pathwayDescription);
            writer.print("</div>\n    <div class=\"show-more\">\n        <a href=\"#\">Show more</a>\n    </div>\n    </div>\n <img src=\"legend.png\" class=\"centeredImage\" alt=\"legend\"> \n</head>\n\n<body>\n\n\n  <div id=\"cy\">\n<script type=\"text/javascript\" src=\"./");
            writer.print(gpmlFileName.replaceFirst("\\.gpml","\\.js"));
            writer.println("\"></script>\n  </div>\n</body>\n</html>");
            writer.close();
            
            //copy legend image
            URL legendURL = getClass().getClassLoader().getResource("medsavant/pathways/cytoscapeJS/displayLegend/legend.png");
            File legendDestination = new File(CACHEFOLDER+OUTPUTDIR+HTMLFOLDER+"legend.png");
            FileUtils.copyURLToFile(legendURL, legendDestination);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
	/**
	 * Build the standard pharmacogenomic conditions to be used when retrieving 
	 * variants for any patient's analysis and store these in a list.
	 * @return a Map of Conditions to be used for all PGx analyses
	 * @throws SQLException
	 */
	private static Map<String, Condition> buildConditionList() throws SQLException {
		Map<String, Condition> queryMap= new HashMap<String, Condition>();
		
                //RUTH - your condition should have DNAid and mutation types.
                //I guess get it working with DNAID for now
                
                
		/* Get all relevant markers for a particular gene and create a
		 * ComboCondition for that set. Then add it to the List. */
			for (String g : PGXDBFunctions.getGenes()) {
				// generate a new query for this gene
				ComboCondition query= new ComboCondition(ComboCondition.Op.OR);
				
				try {
					/* Add all the marker positions for this gene.
					 * NOTE: You can also search for variants using the dbSNP rsID,
					 * however, then you rely on the DB to be up-to-date and annotated
					 * correctly, which is not always the case. It's better to query
					 * variants by chromosomal coordinates. Original code is commented
					 * out below. */
					for (PGXMarker pgxm : PGXDBFunctions.getMarkerInfo(g)) {
						ComboCondition variantCondition= new ComboCondition(ComboCondition.Op.AND);
						variantCondition.addCondition(
							BinaryCondition.equalTo(ts.getDBColumn(BasicVariantColumns.CHROM), pgxm.chromosome));
						variantCondition.addCondition(
							BinaryCondition.equalTo(ts.getDBColumn(BasicVariantColumns.START_POSITION), Integer.parseInt(pgxm.position)));	
						
						query.addCondition(variantCondition);
					}
				} catch (Exception e) {
					e.printStackTrace();
				}
				
				// add this gene-query pair to the list
				queryMap.put(g, query);
			}
		
		return queryMap;		
	}
	
    public int genesInVariants() {
        return this.geneset.size();
    }
    
    
    private void queryVariants(String dnaID, List<String> filterMutationTypes) throws SQLException, RemoteException, SessionExpiredException {
            
        
        /* Take the standard combocondition and AND it to the DNA ID for this
         * individual before submitting for variants. */
        ComboCondition query= new ComboCondition(ComboCondition.Op.AND);
        query.addCondition(
                BinaryCondition.equalTo(ts.getDBColumn(BasicVariantColumns.DNA_ID), dnaID));


        ComboCondition mutationFilter = new ComboCondition(ComboCondition.Op.OR);
        // add conditions for mutation types. Only mutation types in filterMutationTypes will be allowed
        Iterator mutationTypeIterator = filterMutationTypes.iterator();
        String mutationtype;
        while (mutationTypeIterator.hasNext()) {
            mutationtype = (String) mutationTypeIterator.next();
            System.out.println("adding "+mutationtype+" to query");
            mutationFilter.addCondition(
                    //VARIANT_TYPE instead of JANNOVAR EFFECT and then how do you get the variant type out ??
                BinaryCondition.equalTo(ts.getDBColumn(BasicVariantColumns.JANNOVAR_EFFECT.getColumnName()), mutationtype));
        }
        
        query.addCondition(mutationFilter);
        
        /* Once query is built, run it on the remote server. */
        List<Variant> retrievedVariants= runRemoteQuery(query);			
        
        //extract genes from variants for pathway enrichment analysis
        extractGenesFromVariants(retrievedVariants);
        
            
    }
	
    
    private void extractGenesFromVariants(List<Variant> retrievedVariants) {
        this.geneset = new HashSet<String>();
        Iterator<Variant> it = retrievedVariants.iterator();
        String currentGene;
        String[] multipleGenes;
        Variant v;
        while (it.hasNext()) {
            v = it.next(); 
            currentGene = v.getGene();
            multipleGenes = currentGene.split(",");
            //assumption: multiple genes will be in the form "CYP2D6(dist=1145),CYP2D7P1(dist=8186)"
            if (multipleGenes.length > 1) {
                for (int i = 0; i < multipleGenes.length; i++) {
                    currentGene = multipleGenes[i].split("\\(")[0];
                    geneset.add(currentGene);
                    System.out.println("gene: "+currentGene+ " effect: "+v.getMutationSymbols());
                }
            }
            else {
                geneset.add(currentGene);
                System.out.println("gene: "+currentGene+ " effect: "+v.getMutationType());
            }
            
        }
    }
    
    
    
    private List<Variant> runRemoteQuery(Condition query) throws SQLException, RemoteException, SessionExpiredException{
            List<Variant> output= new LinkedList<Variant>();

            /* For each query, a VariantIterator will be returned. When the Iterator
            * is null, stop getting more VariantIterators. Iterate while
            * this object hasNext() and store the Variant objects in a List of
            * Variant objects. Variants are retrieved in chunks based on a request
            * limit offset to allow for a cancellation. */
            Condition[][] conditionMatrix= new Condition[1][1];
            conditionMatrix[0][0]= query;

            int position= 0;
            // initiate VariantIterator for first batch
            List<Object[]> rows= vma.getVariants(LoginController.getInstance().getSessionID(),
                    ProjectController.getInstance().getCurrentProjectID(),
                    ReferenceController.getInstance().getCurrentReferenceID(),
                    conditionMatrix, position, DB_VARIANT_REQUEST_LIMIT);		
            VariantIterator variantIterator= new VariantIterator(rows, ProjectController.getInstance().getCurrentAnnotationFormats());
            while (variantIterator.hasNext()) {
                    // add all the variants to the list from the current batch
                    while (variantIterator != null && variantIterator.hasNext()) {
                            output.add(variantIterator.next());
                    }

                    // increment the request limit
                    position += DB_VARIANT_REQUEST_LIMIT;

                    // Get the next batch 
                    rows= vma.getVariants(LoginController.getInstance().getSessionID(),
                            ProjectController.getInstance().getCurrentProjectID(),
                            ReferenceController.getInstance().getCurrentReferenceID(),
                            conditionMatrix, position, DB_VARIANT_REQUEST_LIMIT);		
                    variantIterator= new VariantIterator(rows, ProjectController.getInstance().getCurrentAnnotationFormats());
            }

            return output;
    }
    
    
    
	
    
    
    
    
    
    
    
    
    
    
    
    
    /**
     * write CSS file for pathway visualization in browser
     * @param fileName 
     */
    private void writeCSS(String fileName) {
        try {
            PrintWriter writer = new PrintWriter(fileName);
            writer.println("body { \n  font: 14px helvetica neue, helvetica, arial, sans-serif;\n}\n\ndiv.text-container {\n\n    text-align: center;\n    margin: 0 auto;\n    width: 75%;    \n\n}\n\n.text-content{\n    line-height: 1em;\n\n    text-align: left;\n\n}\n\n.short-text {\n    overflow: hidden;\n    height: 2em;\n}\n\n.full-text{\n    height: auto;\n}\n\nh1 {\n    font-size: 24px;   \n}\n\n.show-more {\n    padding: 10px 0;\n    text-align: center;\n}\n\n#cy {\n  height: 100%;\n  width: 100%;\n  position: absolute;\n}\n.centeredImage\n{\n    display: block;\n    margin-left: auto;\n    margin-right: auto;\n}\n");
            writer.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    /**
     * Take a node and add its contents to the geneproduct list
     *  or the non geneproduct list. Put it in the position map if
     *  its position is specified. Give it attributes so that it can be
     *  styled correctly
     * @param node
     * @param nonGeneProductNodes
     * @param geneProducts
     * @param positionMap
     * @return 
     */
    private String processNode(Node node, HashMap<String,String> nonGeneProductNodes, HashMap<String,String> geneProducts, HashMap<String,String> positionMap) {
        String nodeTextStart = "      { data: { id: '";
        String nodeTextName = ", name: '";
        String singleQuote = "'";
        String colon = ":";
        String space = " ";
        String nodeTextPlacementNode = ", placementNode: true";
        String isMutatedText = ", isMutated: true";
        String nodeTextWeight = ", weight: ";
        String nodeTextHeight = ", height: ";
        String nodeTextEnd = "} }";
        String startCurlyBrace = "{", endCurlyBrace = "}";
        String x = "x: ", y = ", y: ";
        String width = ", width: ", height = ", height: ";

        
        NamedNodeMap attributes = node.getAttributes();
        String nodeText = "";
        String graphID = "", geneSymbol = "";
        String shapeAttribute = ", shape: ";
        String labelAttribute = ", isLabel: true";
        String inhibitorAttribute = ", isInhibitor: true";
        String metaboliteAttribute = ", isMetabolite: true";
        String parentAttribute = ", parent: '";
        String elbowText=", elbow: true";
        String anchorText = ", anchor: true";
        boolean addToNonGeneProductNodes = false, addToGeneProductNodes = false;
        
        
        //if shape is double lined, call processNode twice
        //once with slightly smaller width/height, once with slightly larger
        if (node.getNodeName().equals("Shape")) {
            if (node.hasChildNodes()) {
                String nullString = "null";
                Node childNode = node.getFirstChild();
                while (childNode!=null && !childNode.getNodeName().equals("Attribute")) {
                    childNode = childNode.getNextSibling();
                }
                if (childNode!=null && childNode.getAttributes().getNamedItem("Key")!=null && childNode.getAttributes().getNamedItem("Key").getNodeValue().equals("org.pathvisio.DoubleLineProperty")) {
                    Node largerCopy = node.cloneNode(true), smallerCopy = node.cloneNode(true);
                    //make unique graphIDs
                    if (node.getAttributes().getNamedItem("GraphId")!=null) {
                        //give smaller node a unique graphID
                        smallerCopy.getAttributes().getNamedItem("GraphId").setNodeValue(smallerCopy.getAttributes().getNamedItem("GraphId").getNodeValue()+"_smaller_outline");
                        //designated larger copy as the smaller copy's parent so that they move together
                        ((Element) smallerCopy).setAttribute("parentID",node.getAttributes().getNamedItem("GraphId").getNodeValue());
                        ((Element) largerCopy).setAttribute("parentID",nullString);
                        
                    }
                    
                    //remove attribute child so that there are not infinite nested double lines
                    this.removeChildAttribute(largerCopy);
                    this.removeChildAttribute(smallerCopy);
                    
                    Node largerGraphics = this.getShapeGraphicsNode(largerCopy);
                    Node smallerGraphics = this.getShapeGraphicsNode(smallerCopy);
                    
                    //compound nodes cannot have custom sizes - remove size attributes of largeCopy
                    ((Element) largerGraphics).removeAttribute("Height");
                    ((Element) largerGraphics).removeAttribute("Width");
                    //change width/height
                    ((Element) smallerGraphics).setAttribute("Height", (Double.parseDouble(((Element) smallerGraphics).getAttribute("Height"))-4.0)+"");
                    ((Element) smallerGraphics).setAttribute("Width", (Double.parseDouble(((Element) smallerGraphics).getAttribute("Width"))-4.0)+"");
                    
                    return this.processNode(largerCopy, nonGeneProductNodes, geneProducts, positionMap) + "\n" + this.processNode(smallerCopy, nonGeneProductNodes, geneProducts, positionMap);
                }
            }
        }
        
        
            
            nodeText+=nodeTextStart;
            
            
            graphID = getGraphRefID(node, nonGeneProductNodes, geneProducts, positionMap).trim();
            nodeText+=graphID+singleQuote;
            String nullString = "null";
            if (attributes.getNamedItem("parentID")==null) {
                String superparentString = "superparent'";
                nodeText+=parentAttribute+superparentString;
            }
            else if (!attributes.getNamedItem("parentID").getNodeValue().equals(nullString)) {
                nodeText+=parentAttribute+attributes.getNamedItem("parentID").getNodeValue()+singleQuote;
            }
            
            nodeText+=nodeTextName;
            
            if (attributes.getNamedItem("TextLabel")!=null) {
                geneSymbol = attributes.getNamedItem("TextLabel").getNodeValue().trim().replace("\n"," ");
                nodeText+=geneSymbol;
            }
            nodeText+= singleQuote;
            if (this.geneset.contains(geneSymbol)) {
                nodeText+=isMutatedText;
            }
            
            if (attributes.getNamedItem("Type")==null || (!attributes.getNamedItem("Type").getNodeValue().equals("GeneProduct") && !attributes.getNamedItem("Type").getNodeValue().equals("Protein") && !attributes.getNamedItem("Type").getNodeValue().equals("Rna"))) {
                addToNonGeneProductNodes = true;
            }
            else if (attributes.getNamedItem("Type").getNodeValue().equals("GeneProduct") || attributes.getNamedItem("Type").getNodeValue().equals("Protein") || attributes.getNamedItem("Type").getNodeValue().equals("Rna")) {
                addToGeneProductNodes = true;
            }
            
            if (attributes.getNamedItem("isPlacementNode")!=null) {
                nodeText+=nodeTextPlacementNode;
            }
            
            if (attributes.getNamedItem("Elbow")!=null) {
                nodeText+=elbowText;
            }
            if (attributes.getNamedItem("Anchor")!=null) {
                nodeText+=anchorText;
            }
            
            if (node.getAttributes().getNamedItem("Type")!=null) {
                if (node.getAttributes().getNamedItem("TextLabel")!=null && node.getAttributes().getNamedItem("TextLabel").getNodeValue().toLowerCase().contains("inhibitor")) {
                    nodeText+=inhibitorAttribute;
                }
                else if (node.getAttributes().getNamedItem("Type").getNodeValue().equals("Metabolite")) {
                    nodeText+=metaboliteAttribute;
                }
                
            }
            
            Node graphics = node.getFirstChild();
            while (graphics!=null && !graphics.getNodeName().equals("Graphics")) {
                graphics = graphics.getNextSibling();
            }
            
            if (graphics!= null ) {
                
                attributes = graphics.getAttributes();


                if (node.getNodeName().equals("Label")) {
                    nodeText+=labelAttribute;
                }
                else if (node.getNodeName().equals("Shape")) {
                    nodeText+=shapeAttribute;
                    if (attributes.getNamedItem("ShapeType")!=null) {
                        String shapetype = attributes.getNamedItem("ShapeType").getNodeValue();
                        if (shapetype.equals("RoundedRectangle")) {
                            nodeText+="'roundrectangle'";
                        }
                        else if (shapetype.equals("Oval")) {
                            nodeText+="'ellipse'";
                        }
                        else {
                            nodeText+="'"+shapetype.toLowerCase().trim()+"'";
                        }
                    }
                    else {
                        nodeText+="'rectangle'";
                    }
                }

                //width/height
                if (attributes.getNamedItem("Width")!=null && attributes.getNamedItem("Height")!=null) {
                    String w = attributes.getNamedItem("Width").getNodeValue();
                    String h = attributes.getNamedItem("Height").getNodeValue();
                    nodeText+=width + w + height + h;
                }
            }
            
            
            nodeText+=nodeTextEnd;
            
            //position
            if (attributes.getNamedItem("CenterX")!=null && attributes.getNamedItem("CenterY")!=null) {
                String xpos = attributes.getNamedItem("CenterX").getNodeValue();
                String ypos = attributes.getNamedItem("CenterY").getNodeValue();
                this.positionMap.put(graphID,singleQuote+graphID+singleQuote+colon+space+startCurlyBrace+x+xpos+y+ypos+endCurlyBrace);
            }
            else if (attributes.getNamedItem("X")!=null && attributes.getNamedItem("Y")!=null) {
                String xpos = attributes.getNamedItem("X").getNodeValue();
                String ypos = attributes.getNamedItem("Y").getNodeValue();
                this.positionMap.put(graphID,singleQuote+graphID+singleQuote+colon+space+startCurlyBrace+x+xpos+y+ypos+endCurlyBrace);
            }
            
            
        
        
        if (addToNonGeneProductNodes) {
            nonGeneProductNodes.put(graphID, nodeText);
        }
        else if (addToGeneProductNodes) {
            geneProducts.put(graphID,nodeText);
        }
        return nodeText;
    }
    /**
     * Retrieve graphics information child node of a shape node.
     * @param node
     * @return 
     */
    private Node getShapeGraphicsNode(Node node) {
        Node childNode = node.getFirstChild();
        while (childNode!=null && !childNode.getNodeName().equals("Graphics")) {
            childNode = childNode.getNextSibling();
        }
        return childNode;
    }
    /**
     * Remove attribute child node of given node
     * @param node 
     */
    private void removeChildAttribute(Node node) {
        Node childNode = node.getFirstChild();
        while (childNode!=null && !childNode.getNodeName().equals("Attribute")) {
            childNode = childNode.getNextSibling();
        }
        node.removeChild(childNode);
    }
    /**
     * Create string representing edge specifications, for inclusion in
     *  javascript file. Extra processing ensures that elbow edges are
     *  correctly stored.
     * @param edgeNode
     * @param nonGeneProductNodes
     * @param geneProducts
     * @return 
     */
    private String processEdge(Node edgeNode, HashMap<String,String> nonGeneProductNodes, HashMap<String,String> geneProducts) {
        boolean arrowHead = true;
        Node target, source, nodeChild = edgeNode.getFirstChild();
        
        //arbitrarily assign points as Source and Target of edge
        while (!nodeChild.getNodeName().equals("Point")) {
            nodeChild = nodeChild.getNextSibling();
        }
        source = nodeChild;
        nodeChild = nodeChild.getNextSibling();
        while (!nodeChild.getNodeName().equals("Point")) {
            nodeChild = nodeChild.getNextSibling();
        }
        target = nodeChild;
        
        
        
        String sourcetext = "      { data: { source: '";
        String targettext = "', target: '";
        String arrowheadtext = "', noArrowHead: true";
        String singlequote = "'", endlinetext = " } },", linebreak = "\n";
        String edgeText="";
        
        //if there are multiple points, no anchors/elbow.
        Node nextPoint;
        nodeChild = nodeChild.getNextSibling();
        while (nodeChild!=null && !nodeChild.getNodeName().equals("Point")) {
            nodeChild = nodeChild.getNextSibling();
        }
        nextPoint = nodeChild;
        if (nextPoint!=null) {
            // determine startelbowdir
            if (edgeNode.getAttributes().getNamedItem("ConnectorType")!=null && edgeNode.getAttributes().getNamedItem("ConnectorType").getNodeValue().equals("Elbow")) {
                boolean startVertical = true;
                Node prevPoint = target;
                //make first source/target edge
                edgeText += getElbowEdgeText(edgeNode,source, target, nonGeneProductNodes, geneProducts, startVertical);
                nodeChild = nextPoint;
                startVertical = !startVertical;
                while (nextPoint!=null) {
                    //make edge from prevpoint to nextpoint
                    edgeText+=linebreak+getElbowEdgeText(edgeNode,prevPoint, nextPoint, nonGeneProductNodes, geneProducts, startVertical);
                    prevPoint = nextPoint;
                    startVertical = !startVertical;
                    nodeChild = nodeChild.getNextSibling();
                    while (nodeChild!=null && !nodeChild.getNodeName().equals("Point")) {
                        nodeChild = nodeChild.getNextSibling();
                    }
                    nextPoint = nodeChild;
                }
                
                
                
                
            }
            //if nto elbow, just generate edges
            else {
                Node prevPoint = target;
                //make first source/target edge
                edgeText+=makeEdgeText(edgeNode,source, target, nonGeneProductNodes, geneProducts);
                nodeChild = nextPoint;
                while (nextPoint!=null) {
                    //make edge from prevpoint to nextpoint
                    edgeText+=linebreak+makeEdgeText(edgeNode,prevPoint, nextPoint, nonGeneProductNodes, geneProducts);
                    prevPoint = nextPoint;
                    nodeChild = nodeChild.getNextSibling();
                    while (nodeChild!=null && !nodeChild.getNodeName().equals("Point")) {
                        nodeChild = nodeChild.getNextSibling();
                    }
                    nextPoint = nodeChild;
                }
            }
            
        }
        else {
            //correct assignment of points as Source and Target of edge - targets have arrowheads
            if (source.getAttributes().getNamedItem("ArrowHead")!=null) {
                Node temp = target;
                target = source;
                source = temp;
            }
            /*Node anchor = edgeNode.getFirstChild();
            while (anchor!=null && !anchor.getNodeName().equals("Anchor")) {
                anchor = anchor.getNextSibling();
            }
            if (anchor!=null) {
                if (anchor.getAttributes().getNamedItem("GraphID")==null && anchor.getAttributes().getNamedItem("GraphRef")==null && anchor.getAttributes().getNamedItem("Type")==null) {
                    ((Element) anchor).removeAttribute("Type");
                    ((Element) anchor).setAttribute("Anchor","true");
                    String anchorGraphID = getGraphRefID(anchor, nonGeneProductNodes, geneProducts, positionMap);
                    ((Element) anchor).setAttribute("GraphId",anchorGraphID);
                }
                edgeText+=makeEdgeText(edgeNode,source, anchor, nonGeneProductNodes, geneProducts);
                edgeText+=linebreak;
                edgeText+=makeEdgeText(edgeNode,anchor, target, nonGeneProductNodes, geneProducts);
            }
            else {*/
                //check if elbow edge (perpendicular)
                if (edgeNode.getAttributes().getNamedItem("ConnectorType")!=null && edgeNode.getAttributes().getNamedItem("ConnectorType").getNodeValue().equals("Elbow") && this.getElbowCoordinates(source, target, geneProducts, nonGeneProductNodes)!=null) {
                    edgeText += getElbowEdgeText(edgeNode, source, target, nonGeneProductNodes, geneProducts, false);
                }
                else {
                    edgeText += makeEdgeText(edgeNode,source, target, nonGeneProductNodes, geneProducts);
                }
           // }
        }
        if (edgeText.contains("[Point: null]")) {
            System.out.println("null node found");
        }
        return edgeText;
    }
    /**
     * get the text for an elbow edge (made up of 2 edges at a right angle)
     * @param edgeNode
     * @param source
     * @param target
     * @param nonGeneProductNodes
     * @param geneProducts
     * @param startVertical
     * @return 
     */
    private String getElbowEdgeText (Node edgeNode, Node source, Node target, HashMap<String,String> nonGeneProductNodes, HashMap<String,String> geneProducts, boolean startVertical) {
        String elbowEdgeText="";
        String linebreak = "\n";
        Node elbow = makeElbow(edgeNode,source, target, nonGeneProductNodes, geneProducts, startVertical);
        elbowEdgeText += makeEdgeText(edgeNode,source, elbow, nonGeneProductNodes, geneProducts);
        elbowEdgeText += linebreak + makeEdgeText(edgeNode,elbow, target, nonGeneProductNodes, geneProducts);
        return elbowEdgeText;
    }
    /**
     * create elbow edge
     * @param edgeNode
     * @param source
     * @param target
     * @param nonGeneProductNodes
     * @param geneProducts
     * @param startVertical
     * @return 
     */
    private Node makeElbow(Node edgeNode, Node source, Node target, HashMap<String,String> nonGeneProductNodes, HashMap<String,String> geneProducts, boolean startVertical) {
        
        //give elbow positioning
        String elbowGraphID = this.getElbowCoordinates(source, target, geneProducts, nonGeneProductNodes, startVertical);
        //need null check
        String[] elbowCoordinates = elbowGraphID.split(" ");
        //make a node for the elbow anchor
        Element elbow = (Element) source.cloneNode(false);
        elbow.setAttribute("X",elbowCoordinates[0]);
        elbow.setAttribute("Y",elbowCoordinates[1]);
        elbow.setAttribute("GraphId",elbowGraphID);
        elbow.removeAttribute("GraphRef");
        elbow.removeAttribute("Type");
        elbow.setAttribute("Elbow","true");
        if (target.getAttributes().getNamedItem("ArrowHead")!=null && (target.getAttributes().getNamedItem("ArrowHead").getNodeValue().equals("mim-inhibition") || target.getAttributes().getNamedItem("ArrowHead").getNodeValue().equals("TBar"))) {
            elbow.setAttribute("ArrowHead","inhibitorElbow");
        }
        elbow.setAttribute("TextLabel","");
        return (Node) elbow;
        
    }
    /**
     * return x and y coordinate of anchor at edge elbow, where the x-coordinate
     *  belongs to the source and the y-coordinate belongs to the target node
     * @param source
     * @param target
     * @param geneProductNodes
     * @param nonGeneProductNodes
     * @return 
     */
    private String getElbowCoordinates(Node source, Node target, HashMap<String,String> geneProductNodes, HashMap<String,String> nonGeneProductNodes) {
        return getElbowCoordinates(source, target, geneProductNodes, nonGeneProductNodes, false);
    }
    /**
     * return x and y coordinate of anchor at edge elbow, taking into account
     *  whether or not the segment connected to the source node should be
     *  vertical or horizontal
     * @param source
     * @param target
     * @param geneProductNodes
     * @param nonGeneProductNodes
     * @param startVertical
     * @return 
     */
    private String getElbowCoordinates(Node source, Node target, HashMap<String,String> geneProductNodes, HashMap<String,String> nonGeneProductNodes, boolean startVertical) {
        String coordinates = "";
        String sourceID, targetID;  
        if (source.getAttributes().getNamedItem("GraphId")!=null) {
            sourceID = source.getAttributes().getNamedItem("GraphId").getNodeValue();
        }
        else if (source.getAttributes().getNamedItem("GraphRef")!=null) {
            sourceID = source.getAttributes().getNamedItem("GraphRef").getNodeValue();
        }
        else {
            sourceID = null;
        }
        if (target.getAttributes().getNamedItem("GraphId")!=null) {
            targetID = target.getAttributes().getNamedItem("GraphId").getNodeValue();
        }
        else if (target.getAttributes().getNamedItem("GraphRef")!=null) {
            targetID = target.getAttributes().getNamedItem("GraphRef").getNodeValue();
        }
        else {
            targetID = null;
        }
        
        String getXPosID, getYPosID;
        Node getXPosNode, getYPosNode;
        if (startVertical) {
            getXPosID = sourceID;
            getXPosNode = source;
            getYPosID = targetID;
            getYPosNode = target;
        }
        else {
            getXPosID = targetID;
            getXPosNode = target;
            getYPosID = sourceID;
            getYPosNode = source;
        }
        
        String nodeText = positionMap.get(getXPosID);
        if (getXPosID!=null && nodeText!=null) {
            //regex extract
            String regEx = "x:([^,]*),";
            Pattern p;
            Matcher m;
            p=Pattern.compile(regEx);
            m = p.matcher(nodeText);
            if (m.find()) {
                coordinates+=m.group(1).trim();
            }
        }
        else if (getXPosNode.getAttributes().getNamedItem("X")!= null) {
            coordinates+=getXPosNode.getAttributes().getNamedItem("X").getNodeValue();
        }
        else if (getXPosNode.getAttributes().getNamedItem("CenterX") != null) {
            coordinates+=getXPosNode.getAttributes().getNamedItem("CenterX").getNodeValue();
        }
        else {
            return null;
        }
        coordinates+=" ";
        nodeText = positionMap.get(getYPosID);
        if (getYPosID!=null && nodeText!=null) {
            //regex extract
            String regEx = "y:([^}]*)}";
            Pattern p;
            Matcher m;
            p=Pattern.compile(regEx);
            m = p.matcher(nodeText);
            if (m.find()) {
                coordinates+=m.group(1).trim();
            }
        }
        else if (getYPosNode.getAttributes().getNamedItem("Y")!= null) {
            coordinates+=getYPosNode.getAttributes().getNamedItem("Y").getNodeValue();
        }
        else if (getYPosNode.getAttributes().getNamedItem("CenterY")!=null) {
            coordinates+=getYPosNode.getAttributes().getNamedItem("CenterY").getNodeValue();
        }
        else {
            return null;
        }
        return coordinates;
    }
    /**
     * make text representing edge for javascript file
     * @param graphics
     * @param source
     * @param target
     * @param nonGeneProductNodes
     * @param geneProducts
     * @return 
     */
    private String makeEdgeText(Node graphics, Node source, Node target, HashMap<String,String> nonGeneProductNodes, HashMap<String,String> geneProducts) {
        String sourcetext = "      { data: { source: '";
        String targettext = ", target: '";
        String labelText = ", label: '";
        String noArrowHeadText = ", noArrowHead: true";
        String circularArrowHeadText = ", circularArrowHead: true";
        String veeArrowHeadText = ", veeArrowHead: true";
        String triangleArrowHeadText = ", triangleArrowHead: true";
        String teeArrowHeadText = ", teeArrowHead: true";
        String singlequote = "'", endlinetext = " } },";
        String inhibitionString = "mim-inhibition";
        String tbar = "TBar";
        String catalysisString = "mim-catalysis";
        String bindingString = "mim-binding";
        String inhibitorElbowText = "inhibitorElbow";
        String haystackText = ", haystackStyle: true";
        String redColorText = ", color: 'red'";
        
        String edgeText = sourcetext;
        edgeText+= this.getEdgeGraphRefID(source, nonGeneProductNodes, geneProducts, positionMap)+singlequote;
        edgeText+=targettext;
        edgeText+= this.getEdgeGraphRefID(target, nonGeneProductNodes, geneProducts, positionMap)+singlequote;
        
        //MIM arrowhead conventions can be found here: http://discover.nci.nih.gov/mim/mapDesc.html
        
        if (target.getAttributes().getNamedItem("ArrowHead")==null) { //no arrowhead
            edgeText+=noArrowHeadText;
        }
        else {
            edgeText+=labelText;
            
            String arrowHeadType = target.getAttributes().getNamedItem("ArrowHead").getNodeValue();
            String[] arrowHeadLabel = arrowHeadType.split("-");
            edgeText+=arrowHeadLabel[arrowHeadLabel.length-1]+singlequote;
            
            if (arrowHeadType.equals(catalysisString)) {
                edgeText+=circularArrowHeadText;
            }
            else if (arrowHeadType.equals(inhibitionString) || arrowHeadType.equals(tbar)) {
                edgeText+=teeArrowHeadText+redColorText;
            }
            else if (arrowHeadType.equals(bindingString)) {
                edgeText+=veeArrowHeadText;
            }
            else if (arrowHeadType.equals(inhibitorElbowText)) {
                edgeText+=redColorText+noArrowHeadText;
            }
            else {
                edgeText+=triangleArrowHeadText;
            }
        }
        if (graphics.getAttributes().getNamedItem("ConnectorType")!=null && graphics.getAttributes().getNamedItem("ConnectorType").getNodeValue().equals("Elbow")) {
            edgeText+=haystackText;
        }
        edgeText+=endlinetext;
        return edgeText;
    }
    /**
     * Retrieve graphID of node. If node doesn't have one, use positional
     *  information (x_position space y_position) as the ID
     * @param node
     * @param nonGeneProductNodes
     * @param geneProducts
     * @param positionMap
     * @return 
     */
    private String getGraphRefID(Node node, HashMap<String,String> nonGeneProductNodes, HashMap<String,String> geneProducts, HashMap<String,String> positionMap) {
        String graphID;
        boolean newNonGeneProduct = false;
        if (node.getAttributes().getNamedItem("GraphRef")!=null) {
            graphID = node.getAttributes().getNamedItem("GraphRef").getNodeValue().trim();
            
        }
        else if (node.getAttributes().getNamedItem("GraphId")!=null) {
            graphID = node.getAttributes().getNamedItem("GraphId").getTextContent().trim();
        }
        else {
            graphID = getPositionGraphID(node);
        }
        //check centerx and figure out some defautl... maybe position
        if (!geneProducts.containsKey(graphID) && !nonGeneProductNodes.containsKey(graphID)) {
            //check if there is positional information
            if (node.getAttributes().getNamedItem("X")!=null){
                positionMap.put(graphID,"'"+graphID+"': {x: "+node.getAttributes().getNamedItem("X").getNodeValue()+", y: "+node.getAttributes().getNamedItem("Y").getNodeValue()+"}");
            }
            else if (node.getAttributes().getNamedItem("CenterX")!=null) {
                positionMap.put(graphID,"'"+graphID+"': {x: "+node.getAttributes().getNamedItem("CenterX").getNodeValue()+", y: "+node.getAttributes().getNamedItem("CenterY").getNodeValue()+"}");
            }
        }
        return graphID;
    }
    /**
     * Returns string representing graph ID based on position, composed of
     *  x_position space y_position
     * @param node
     * @return 
     */
    private String getPositionGraphID(Node node) {
        Node graphics = node.getFirstChild();
        String graphID = "";
        while (graphics!=null && graphics.getNodeName()!="Graphics"){
            graphics = graphics.getNextSibling();
        }
        if (graphics==null) {
            if (node.getAttributes().getNamedItem("X")!=null){
                graphID = node.getAttributes().getNamedItem("X").getNodeValue() + " " + node.getAttributes().getNamedItem("Y").getNodeValue().trim();
            }
            else if (node.getAttributes().getNamedItem("CenterX")!=null) {
                graphID = node.getAttributes().getNamedItem("CenterX").getNodeValue() + " " + node.getAttributes().getNamedItem("CenterY").getNodeValue().trim();
            }
            else {
                System.out.println("there is no way of assigning an ID to this node.");
            }
        }
        else {
            if (graphics.getAttributes().getNamedItem("X")!=null){
                graphID = graphics.getAttributes().getNamedItem("X").getNodeValue() + " " + graphics.getAttributes().getNamedItem("Y").getNodeValue().trim();
            }
            else if (graphics.getAttributes().getNamedItem("CenterX")!=null) {
                graphID = (graphics.getAttributes().getNamedItem("CenterX").getNodeValue() + " " + graphics.getAttributes().getNamedItem("CenterY").getNodeValue()).trim();
            }
            else {
                System.out.println("there is no way of assigning an ID to this node, but it has graphics.");
            }
        }
        return graphID;
    }
    /**
     * get graph ID during processing of edges (new node is made if the target
     *  or source is not already a node)
     * @param node
     * @param nonGeneProductNodes
     * @param geneProducts
     * @param positionMap
     * @return 
     */
    private String getEdgeGraphRefID(Node node, HashMap<String,String> nonGeneProductNodes, HashMap<String,String> geneProducts, HashMap<String,String> positionMap) {
        if (node == null) {
            System.out.println("tried to find graphID of null node");
        }
        String graphID = getGraphRefID(node,nonGeneProductNodes, geneProducts, positionMap);
        if (!geneProducts.containsKey(graphID) && !nonGeneProductNodes.containsKey(graphID)) {
            ((Element) node).setAttribute("isPlacementNode","true");
            processNode(node,nonGeneProductNodes, geneProducts, this.positionMap);
        }
        else if (!this.positionMap.containsKey(graphID)) {
            String singleQuote = "'";
            String colon = ":";
            String space = " ";
            String startCurlyBrace = "{";
            String endCurlyBrace = "}";
            String x = "x: ", y = ", y: ";
            NamedNodeMap attributes = node.getAttributes();
            if (attributes.getNamedItem("CenterX")!=null && attributes.getNamedItem("CenterY")!=null) {
                String xpos = attributes.getNamedItem("CenterX").getNodeValue();
                String ypos = attributes.getNamedItem("CenterY").getNodeValue();
                this.positionMap.put(graphID,singleQuote+graphID+singleQuote+colon+space+startCurlyBrace+x+xpos+y+ypos+endCurlyBrace);
            }
            else if (attributes.getNamedItem("X")!=null && attributes.getNamedItem("Y")!=null) {
                String xpos = attributes.getNamedItem("X").getNodeValue();
                String ypos = attributes.getNamedItem("Y").getNodeValue();
                this.positionMap.put(graphID,singleQuote+graphID+singleQuote+colon+space+startCurlyBrace+x+xpos+y+ypos+endCurlyBrace);
            }
        }
        return graphID;
    }
    /**
     * extract genes and effects from variant call format file
     * @param filename
     * @param GENE_REGEX
     * @param EFFECT_REGEX
     * @param GENES_NOT_IN_GENESETS_FILE 
     */
    public void readVCF(String filename) {
        //dgenes = new ArrayList<String>();
        //deffects = new ArrayList<String>();
        geneset = new HashSet<String>();
        effectset = new HashSet<String> ();
        Pattern p;
        Matcher m;
        try {
            BufferedReader br = new BufferedReader(new FileReader(filename));
            String line = br.readLine();
            String temp;
            while (line != null) {
                //skip header/preamble lines starting with #
                line = line.trim();
                if (line.charAt(0)!='#') {
                    //8th column has extra information separated by ;
                    p=Pattern.compile(GENE_REGEX);
                    m = p.matcher(line);
                    if (m.find()) {
                        geneset.add(m.group(1));
                        p=Pattern.compile(EFFECT_REGEX);
                        m = p.matcher(line);
                        m.find();
                        effectset.add(m.group(1));
                    }
                }
                line = br.readLine();
            }
            HashSet<String> genesNotInGeneSets = new HashSet<String>(geneset);
            genesNotInGeneSets.removeAll(allgenes);
            PrintWriter writer = new PrintWriter(GENES_NOT_IN_GENESETS_FILE, "UTF-8");
            writer.println("Gene symbol");
            Iterator<String> nonGenesetGenes = genesNotInGeneSets.iterator();
            while (nonGenesetGenes.hasNext()) {
                writer.println(nonGenesetGenes.next());
            }
            writer.close();
            geneset.retainAll(allgenes); //geneset is the set of mutated genes
            
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    /**
     * extract genes and effects from variant call format file
     * @param filename
     * @param GENE_REGEX
     * @param EFFECT_REGEX
     * @param GENES_NOT_IN_GENESETS_FILE 
     */
    public void readVCF(File filename) {
        geneset = new HashSet<String>();
        effectset = new HashSet<String> ();
        Pattern p;
        Matcher m;
        try {
            System.out.println("reached readVCF method");
            BufferedReader br = new BufferedReader(new FileReader(filename));
            System.out.println("able to parse filie into buffered reader");
            String line = br.readLine();
            while (line != null) {
                //skip header/preamble lines starting with #
                line = line.trim();
                if (line.charAt(0)!='#') {
                    //8th column has extra information separated by ;
                    p=Pattern.compile(GENE_REGEX);
                    m = p.matcher(line);
                    if (m.find()) {
                        geneset.add(m.group(1));
                        p=Pattern.compile(EFFECT_REGEX);
                        m = p.matcher(line);
                        m.find();
                        effectset.add(m.group(1));
                    }
                }
                line = br.readLine();
            }
            HashSet<String> genesNotInGeneSets = new HashSet<String>(geneset);
            genesNotInGeneSets.removeAll(allgenes);
            PrintWriter writer = new PrintWriter(GENES_NOT_IN_GENESETS_FILE, "UTF-8");
            writer.println("Gene symbol");
            Iterator<String> nonGenesetGenes = genesNotInGeneSets.iterator();
            while (nonGenesetGenes.hasNext()) {
                writer.println(nonGenesetGenes.next());
            }
            writer.close();
            geneset.retainAll(allgenes); //geneset is the set of mutated genes
            
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    /**
     * perform analysis that produces a p-value to show how likely it is that
     *  a pathway is enriched for genetic variants
     */
    public void hypergeometricTest() {
        //for each pathway
        HashSet<String> commonGenes;
        String pathwayname;
        Iterator pathwaynames = genesets.keySet().iterator();
        //marked items is the number of genes in VCF file
        //samplesize is number of genes in VCF file AND pathway in question
        //populationsize is all genes in genesets
        int markeditems = geneset.size(), samplesize, populationsize = allgenes.size();
        Hypergeometric h;
        testedPathways = new ArrayList<TestedPathway>();
        double p;
        while (pathwaynames.hasNext()) {
            pathwayname = (String) pathwaynames.next();
            //get pathway
            commonGenes = new HashSet<String>((ArrayList<String>) genesets.get(pathwayname));
            samplesize = commonGenes.size();
            commonGenes.retainAll(geneset);
            //see if there are common genes
            if (commonGenes.size() > 0) {
                try {
                    h = new Hypergeometric(samplesize, populationsize, markeditems);
                    p = h.cdf((double) commonGenes.size());
                    //store p-value from hypergeometric test
                    //if p is less than some threshhold
                    //color all genes involved in the pathway
                    testedPathways.add(new TestedPathway(p, pathwayname,commonGenes, pathwayGpmls.get(pathwayname), ((ArrayList<String>) genesets.get(pathwayname)).size()));
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
        
        //comparators.sort by cdf
        Collections.sort(testedPathways, new PathwayEnrichmentProbabilityComparator());
        //comparators sort by pathway name (HUMAN READABLE)
        //Collections.sort(testedpathways, new PathwayNameComparator());
    }
    
    
    /**
     * perform analysis that produces a p-value to show how likely it is that
     *  a pathway is enriched for genetic variants
     * @param WIKIPATHWAYSFOLDER 
     */
    public void hypergeometricWikiPathwaysTest() {
        //for each pathway
        HashSet<String> commonGenes;
        String pathwayname;
        Iterator pathwaynames = genesets.keySet().iterator();
        int samplesize = geneset.size(), markeditems, populationsize = allgenes.size();
        Hypergeometric h;
        testedPathways = new ArrayList<TestedPathway>();
        double p;
        int numPathwayGenes;
        while (pathwaynames.hasNext()) {
            pathwayname = (String) pathwaynames.next();
            numPathwayGenes =((HashMap) genesets.get(pathwayname)).size();
            if (numPathwayGenes >= minPathwayGenesFilter && numPathwayGenes <= maxPathwayGenesFilter) {
                //get pathway
                commonGenes = new HashSet<String>( ( (HashMap<String,Integer>) genesets.get(pathwayname) ).keySet());
                markeditems = commonGenes.size();
                commonGenes.retainAll(geneset);
                //see if there are common genes
                //if (commonGenes.size() > 0) {
                    try {
                        System.out.println("HYPERGEOMETRIC TEST: sample size - "+samplesize+", population size - "+populationsize+", marked items - "+markeditems);
                        h = new Hypergeometric(samplesize, populationsize, markeditems);
                        p = 1.0-h.cdf((double) commonGenes.size());
                        testedPathways.add(new TestedPathway(p, pathwayname,commonGenes, pathwayGpmls.get(pathwayname), numPathwayGenes));
                    }
                    catch (Exception e) {
                        e.printStackTrace();
                    }
                //}
            }
        }
        
        //comparators.sort by cdf
        Collections.sort(testedPathways, new PathwayEnrichmentProbabilityComparator());
        //comparators sort by pathway name (HUMAN READABLE)
        //Collections.sort(testedpathways, new PathwayNameComparator());
    }
    
    public void generateGPMLFile(String pathwayTitle, String pathwayFileName) {
       HashSet<String> commonGenes = new HashSet<String>( ( (HashMap<String,Integer>) genesets.get(pathwayTitle) ).keySet());
       //set geneset - give each pathway geneset i guess beforehand
        commonGenes.retainAll(geneset);
        this.markPathwayGenesInGPML(commonGenes, pathwayFileName, CACHEFOLDER+OUTPUTDIR+GPMLFOLDER, PVALUE_CUTOFF);
                    
    }
    
    /**
     * write list of pathways and p-values to table
     * @param filename 
     */
    public void outputEnrichedGeneList(String filename) {
        System.out.println("Writing pathways and p-values to "+filename);
        try {
            PrintWriter writer = new PrintWriter(filename, "UTF-8");
            writer.println("Gene symbol\tP-value\tGenes");
            TestedPathway t;
            HashSet<String> pathwayGenes;
            Iterator it;
            int size = testedPathways.size();
            for (int i = 0; i < size; i++) {
                t = testedPathways.get(i);
                writer.print(t.getName() + "\t" + t.getP());
                pathwayGenes=t.getGenes();
                it = pathwayGenes.iterator();
                while (it.hasNext()) {
                    writer.print("\t"+it.next());
                }
                writer.println();
            }
            writer.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    /**
     * create genesets from wikipathways GPML files
     * @param pathwayfolder
     * @param OUTFILE 
     */
    /*public void wikipathways2GMT() {
        
        try {
            Document doc;
            XPath xpath;
            String xpathexpression;
            NodeList nodes;
            String geneNodeName;
            String[] geneNames;
            PrintWriter writer = new PrintWriter(new File(this.getClass().getResource("geneSet/wikipathways.gmt").getPath()));
            PathMatchingResourcePatternResolver resolver = new PathMatchingResourcePatternResolver();
            Resource[] gpmlFileResources = resolver.getResources("classpath*:**//*}}{{*.gpml");
            //get number of folders in fileEntry list
            int numPathways = gpmlFileResources.length;
            this.pathwayTitles = new String[numPathways];
            this.pathwayHtmlFileNames = new String[numPathways];
            this.pathwayGpmlFileNames = new String[numPathways];
            int pathwayIndex = 0;
            
            File fileEntry;
            for (int i = 0; i < numPathways; i++) {
                fileEntry = gpmlFileResources[i].getFile();
                pathwayGpmlFileNames[pathwayIndex] = fileEntry.getName();
                pathwayHtmlFileNames[pathwayIndex] = OUTPUTDIR+HTMLFOLDER+fileEntry.getName().replaceFirst("\\.gpml","\\.html");
                doc = DocumentBuilderFactory.newInstance().newDocumentBuilder().parse(new InputSource(fileEntry.getPath()));
                xpath = XPathFactory.newInstance().newXPath();
                
                writer.print(fileEntry.getName());
                
                xpathexpression = "/Pathway/@Name";
                nodes = (NodeList) xpath.evaluate(xpathexpression,doc,XPathConstants.NODESET);
                if (nodes.getLength() >0 && nodes.item(0).getNodeType()==Node.ATTRIBUTE_NODE){
                    writer.print("\t"+((Attr) nodes.item(0)).getValue());
                    pathwayTitles[pathwayIndex] = ((Attr) nodes.item(0)).getValue();
                }
                else {
                    System.out.println("unable to find pathway name in "+fileEntry.getName());
                    pathwayTitles[pathwayIndex] = "No pathway title";
                }
                pathwayIndex++;
                xpathexpression = "/Pathway/Comment[@Source='WikiPathways-description']";
                nodes = (NodeList) xpath.evaluate(xpathexpression,doc,XPathConstants.NODESET);
                if (nodes.getLength() >0){
                    writer.print("|"+nodes.item(0).getTextContent().replaceAll("\\n"," "));
                }
                HashSet<String> singlePathwayGenes = new HashSet<String>();
                
                xpathexpression = "/Pathway/DataNode[@Type='GeneProduct']";
                nodes = (NodeList) xpath.evaluate(xpathexpression,doc,XPathConstants.NODESET);
                for (int counter = 0; counter < nodes.getLength(); counter++) {
                    if (nodes.item(counter).getNodeType()==Node.ELEMENT_NODE) {
                        geneNodeName = ((Element) nodes.item(counter)).getAttribute("TextLabel");
                        geneNames = geneNodeName.split("[\\/]|\\s");
                        for (int counter2 = 0; counter2<geneNames.length;counter2++) {
                            geneNames[counter2] = geneNames[counter2].trim();
                            //make sure there are no spaces and no lower case letters (indicative of a non-gene symbol)
                            Pattern p = Pattern.compile("[^A-Z0-9-]");
                            Matcher m = p.matcher(geneNames[counter2]);
                            if (geneNames[counter2].length() != 0 && !m.find()) {
                                singlePathwayGenes.add(geneNames[counter2]);
                            }
                        }
                    }
                }
                xpathexpression = "/Pathway/DataNode[@Type='Protein']";
                nodes = (NodeList) xpath.evaluate(xpathexpression,doc,XPathConstants.NODESET);
                for (int counter = 0; counter < nodes.getLength(); counter++) {
                    if (nodes.item(counter).getNodeType()==Node.ELEMENT_NODE) {
                        geneNodeName = ((Element) nodes.item(counter)).getAttribute("TextLabel");
                        geneNames = geneNodeName.split("[\\/]|\\s");
                        for (int counter2 = 0; counter2<geneNames.length;counter2++) {
                            geneNames[counter2] = geneNames[counter2].trim();
                            //make sure there are no spaces and no lower case letters (indicative of a non-gene symbol)
                            Pattern p = Pattern.compile("[^A-Z0-9-]");
                            Matcher m = p.matcher(geneNames[counter2]);
                            if (geneNames[counter2].length() != 0 && !m.find()) {
                                singlePathwayGenes.add(geneNames[counter2]);
                            }
                        }
                    }
                }
                xpathexpression = "/Pathway/DataNode[@Type='Rna']";
                nodes = (NodeList) xpath.evaluate(xpathexpression,doc,XPathConstants.NODESET);
                for (int counter = 0; counter < nodes.getLength(); counter++) {
                    if (nodes.item(counter).getNodeType()==Node.ELEMENT_NODE) {
                        geneNodeName = ((Element) nodes.item(counter)).getAttribute("TextLabel");
                        geneNames = geneNodeName.split("[\\/]|\\s");
                        for (int counter2 = 0; counter2<geneNames.length;counter2++) {
                            geneNames[counter2] = geneNames[counter2].trim();
                            //make sure there are no spaces and no lower case letters (indicative of a non-gene symbol)
                            Pattern p = Pattern.compile("[^A-Z0-9-]");
                            Matcher m = p.matcher(geneNames[counter2]);
                            if (geneNames[counter2].length() != 0 && !m.find()) {
                                singlePathwayGenes.add(geneNames[counter2]);
                            }
                        }
                    }
                }
                
                //print genes into genesets file
                Iterator it = singlePathwayGenes.iterator();
                while (it.hasNext()) {
                    writer.print("\t"+it.next());
                }
                
                writer.println();
            }
            writer.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }*/
    
    
    public ImageIcon getPathwayImage(String pathwayName) {
        return this.pathwayLinks.get(pathwayName);
    }
    
    
    /**
     * increase border thickness of nodes in GPML files that have variants
     *  associated with their gene product
     * @param WIKIPATHWAYSFOLDER
     * @param commonGenes
     * @param pathwayName
     * @param GPMLOUTPUTFOLDER
     * @param PVALUE_CUTOFF 
     */
    public void markPathwayGenesInGPML (HashSet<String> commonGenes, String pathwayFileName, String GPMLOUTPUTFOLDER, double PVALUE_CUTOFF) {
        String gene,geneMatch;
        try {
            //ClassLoader classloader = this.getClass().getClassLoader();
            //PathMatchingResourcePatternResolver resolver = new PathMatchingResourcePatternResolver(classloader);
            //Resource[] gpmlFile = resolver.getResources("classpath*:medsavant/pathways/wikipathwaysGPML/"+pathwayFileName);
            //Document doc = DocumentBuilderFactory.newInstance().newDocumentBuilder().parse(new InputSource(gpmlFile[0].getFile().getPath()));
            System.out.println(this.getClass().getClassLoader().getResourceAsStream("medsavant/pathways/wikipathwaysGPML"));
            System.out.println("pathway filename is: "+pathwayFileName);
            Document doc = DocumentBuilderFactory.newInstance().newDocumentBuilder().parse(new InputSource(new BufferedReader(new InputStreamReader(this.getClass().getClassLoader().getResourceAsStream(WIKIPATHWAYSFOLDER+pathwayFileName)))));
            
            XPath xpath = XPathFactory.newInstance().newXPath();
            
            Iterator enrichedGenes = commonGenes.iterator();
            
            //if there aren't any enriched genes, copy GPML file over as-is
            if (!enrichedGenes.hasNext()) {
                URI uri = this.getClass().getClassLoader().getResource(WIKIPATHWAYSFOLDER+pathwayFileName).toURI();
                final Map<String, String> env = new HashMap<>();
                final String[] array = uri.toString().split("!");
                final FileSystem fs = FileSystems.newFileSystem(URI.create(array[0]), env);
                System.out.println("uri: "+uri);
                 Path path = java.nio.file.Files.copy( 
                       fs.getPath(array[1]), 
                       new java.io.File(CACHEFOLDER+OUTPUTDIR+GPMLFOLDER+pathwayFileName).toPath(),
                       java.nio.file.StandardCopyOption.REPLACE_EXISTING,
                       java.nio.file.StandardCopyOption.COPY_ATTRIBUTES,
                       java.nio.file.LinkOption.NOFOLLOW_LINKS );
                 
                fs.close();
            }
            
            while (enrichedGenes.hasNext()) {
                
                gene = (String) enrichedGenes.next();
                gene = gene.trim();
                String xpathexpression = "/Pathway/DataNode[contains(@TextLabel,'"+gene+"') and (@Type='GeneProduct' or @Type='Protein' or @Type='Rna')]";
                NodeList nodes = (NodeList) xpath.evaluate(xpathexpression,doc,XPathConstants.NODESET);
                boolean found = false;
                for (int counter = 0; counter < nodes.getLength(); counter++) {
                    if (nodes.item(counter).getNodeType()==Node.ELEMENT_NODE) {
                        //check what textlabel is (regex)
                        geneMatch = ((Element) nodes.item(counter)).getAttribute("TextLabel");
                        //get graphic element
                        
                        Pattern p = Pattern.compile("(^|\\s+|[/ '\"]+)"+gene+"(\\s+|[/ '\"]+|$)");
                        Matcher m = p.matcher(geneMatch);
                        if (m.find()) {
                            NodeList graphicComponents = ((Element) nodes.item(counter)).getElementsByTagName("Graphics");
                            ((Element) graphicComponents.item(0)).setAttribute("LineThickness","5.0");
                            found = true;
                        }
                        
                    }
                }
                if (!found) {
                    System.out.println("could not find "+gene+" in "+pathwayFileName + ", nodelength: "+nodes.getLength());
                    //also not fidning things with slashes
                }
            }
            Transformer xformer = TransformerFactory.newInstance().newTransformer();
            xformer.transform(new DOMSource(doc), new StreamResult(new File(CACHEFOLDER+OUTPUTDIR+GPMLFOLDER+pathwayFileName)));
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}
/**
 * allows testedpathways to be sorted by p-value
 * @author ruthgrace
 */
class PathwayEnrichmentProbabilityComparator implements Comparator<TestedPathway> {
    @Override
    public int compare(TestedPathway a, TestedPathway b) {
        if (a.getP() == b.getP()) {
            return 0;
        }
        if (a.getP() > b.getP()) {
            return 1;
        }
        return -1;
    }
    
}
/**
 * allows tested pathways to be sorted by pathway title
 * @author ruthgrace
 */
class PathwayNameComparator implements Comparator<TestedPathway> {
    private final String NAME_REGEX = "\t([^\t])";
    @Override
    public int compare(TestedPathway a, TestedPathway b) {
        String readableNameA,readableNameB;
        Pattern p = Pattern.compile(NAME_REGEX);
        Matcher m = p.matcher(a.getName());
        m.find();
        readableNameA = m.group(1);
        m = p.matcher(b.getName());
        m.find();
        readableNameB = m.group(1);
        return readableNameA.compareToIgnoreCase(readableNameB);
    }
}