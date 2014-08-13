/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package medsavant.pathways;


/**
 * This class represents a pathway that has undergone hypergeometric testing.
 * It has a p-value, a pathway name, and a set of genes that belong to the pathway.
 * @author ruthgrace
 */

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;

public class TestedPathway {
    //p-value for this pathway
    private double p;
    //pathway name
    private String name;
    //gene symbols of genes both in pathway and in the sample
    private HashSet<String> genes;
    //file names associated with this pathway
    private String htmlFilePath;
    private String pngFilePath;
    private String gpmlFilePath;
    //total number of genes in pathway
    private int totalGenes;
    //indices of columns in table
    public static int PATHWAYNAMEINDEX = 1;
    public static int PVALUEINDEX = 0;
    public static int GENESINDEX = 2;
    public static int NUMGENESINDEX = 3;
    //q value for this pathway (FDR)
    private double q;
    
    /**
     * Initialize the tested pathway object, which represents the result of
     *  testing a pathway by the hypergeometric test.
     * @param p p-value
     * @param name pathway title
     * @param genes gene symbols of genes associated with variants in the
     *  sample, which appear in this pathway.
     * @param gpmlPath name of associated GPML file
     * @param totalGenes total number of genes in this pathway
     */
    public TestedPathway(double p, String name, HashSet<String> genes, String gpmlPath, int totalGenes) {
        this.p = p;
        this.name = name;
        this.genes = new HashSet<String>(genes);
        this.gpmlFilePath = gpmlPath;
        this.htmlFilePath = gpmlFilePath.replaceFirst("\\.gpml","\\.html");
        this.pngFilePath = gpmlFilePath.replaceFirst("\\.gpml","\\.png");
        this.totalGenes = totalGenes;
        this.q = 1.0;
    }
    
    /**
     * Retrieve p-value for this pathway
     * @return double p-value
     */
    public double getP() {
        return this.p;
    }
    
    /**
     * Retrieve pathway title
     * @return String pathway title
     */
    public String getName() {
        return this.name;
    }
    
    /**
     * Retrieve genes for this pathway
     * @return hashset of string genes
     */
    public HashSet<String> getGenes() {
        return this.genes;
    }
    
    /**
     * Return an object list for this tested pathway, with all the information
     *  necessary to put this pathway in the results table
     * @return object array of p-value, name, genes, and total genes
     */
    public Object[] getObjectList() {
        if (genes.size() == 0) {
            return null;
        }
        Object[] objectArray = new Object[4];
        objectArray[0] = this.getP();
        objectArray[1] = this.name;
        objectArray[2] = genesToString();
        objectArray[3] = this.totalGenes;
        return objectArray;
    }
    
    /**
     * Return an object list for this tested pathway, with all the information
     *  necessary to put this pathway in the results table, using the
     *  Bonferroni correction.
     * @param numTests number of pathways tested
     * @return object array of corrected p-value, pathway name, genes, and
     *  total genes
     */
    public Object[] getObjectList(int numTests) {
        if (genes.size() == 0) {
            return null;
        }
        Object[] objectArray = new Object[4];
        objectArray[0] = this.getCorrectedPValue(numTests);
        objectArray[1] = this.name;
        objectArray[2] = genesToString();
        objectArray[3] = this.totalGenes;
        return objectArray;
    }
    
    /**
     * Convert a list of tested pathways into object arrays to be passed to the
     *  results table (with Bonferroni multiple test correction)
     * @param list of tested pathways
     * @return list of object with table data
     */
    public static List<Object[]> convertToObjectListBonferroni(List<TestedPathway> list) {
        List<Object[]> newList = new ArrayList<Object[]>();
        Iterator it = list.iterator();
        Object[] objectArray;
        int numTests = list.size();
        while (it.hasNext()) {
            objectArray = ((TestedPathway) it.next()).getObjectList(numTests);
            if (objectArray!=null) {
                newList.add(objectArray);
                
            }
        }
        return newList;
    }
    
    /**
     * Convert a list of tested pathways into object arrays to be passed to the
     *  results table (with no multiple test correction -- never used)
     * @param list of tested pathways
     * @return list of object with table data
     */
    public static List<Object[]> convertToObjectListNoCorrection(List<TestedPathway> list) {
        List<Object[]> newList = new ArrayList<Object[]>();
        Iterator it = list.iterator();
        Object[] objectArray;
        int numTests = list.size();
        while (it.hasNext()) {
            objectArray = ((TestedPathway) it.next()).getObjectList(numTests);
            if (objectArray!=null) {
                newList.add(objectArray);
                
            }
        }
        return newList;
    }
    
    /**
     * Calculate Q-values of all tested pathways for Benjamini Hochberg
     *  multiple test correction and assign it to the tested pathways
     * @param list of tested pathways
     */
    private static void assignQValues(List<TestedPathway> list) {
        Collections.sort(list, new PathwayEnrichmentProbabilityComparator());
        int size = list.size();
        for (int i = 0; i < size; i++) {
            list.get(i).setQ(((double) size)*list.get(i).getP()/((double) i));
        }
    }
    
    /**
     * Return an object list for this tested pathway, with all the information
     *  necessary to put this pathway in the results table, using the
     *  Benjamini-Hochberg correction.
     * @param fdrCutoff the FDR cutoff for the multiple test correction
     * @return object array of Q-value, pathway name, genes, and
     *  total genes
     */
    private Object[] getBHObjectList(double fdrCutoff) {
        if (genes.size() == 0) {
            return null;
        }
        Object[] objectArray = new Object[4];
        objectArray[0] = this.getQ();
        if ((double) objectArray[0] > fdrCutoff) {
            return null;
        }
        objectArray[1] = this.name;
        objectArray[2] = genesToString();
        objectArray[3] = this.totalGenes;
        return objectArray;
    }
    
    /**
     * Convert a list of tested pathways into object arrays to be passed to the
     *  results table (with Benjamini-Hochberg correction)
     * @param list of tested pathways
     * @param fdrCutoff FDR cutoff for multiple test correction
     * @return list of object with table data
     */
    public static List<Object[]> convertToObjectListBH(List<TestedPathway> list, double fdrCutoff) {
        // MAKE BENJAMINI HOCHBERG P-VALUES
        assignQValues(list);
        List<Object[]> newList = new ArrayList<Object[]>();
        Iterator it = list.iterator();
        Object[] objectArray;
        int numTests = list.size();
        while (it.hasNext()) {
            objectArray = ((TestedPathway) it.next()).getBHObjectList(fdrCutoff);
            if (objectArray!=null) {
                newList.add(objectArray);
                
            }
        }
        return newList;
    }
    
    /**
     * Consolidate information about all the tested pathways into a two
     *  dimensional array of String
     * @param list of tested pathways
     * @return String array of pathway titles, html file paths, gpml file paths,
     *  and p-values.
     */
    public static String[][] getPathwayInfoArray (List<TestedPathway> list) {
        Iterator it = list.iterator();
        ArrayList<String> pathwayTitles = new ArrayList<String>();
        ArrayList<String> pathwayHtmlFileNames = new ArrayList<String>();
        ArrayList<String> pathwayGpmlFileNames = new ArrayList<String>();
        ArrayList<String> pathwayPValues = new ArrayList<String>();
        TestedPathway pathway;
        
        while (it.hasNext()) {
            pathway = (TestedPathway) it.next();
            if (!pathway.noGenes()) {
                pathwayTitles.add(pathway.getTitle());
                pathwayHtmlFileNames.add(pathway.getHtml());
                pathwayGpmlFileNames.add(pathway.getGpml());
                pathwayPValues.add(String.format("%1.2e",pathway.getP()));
            }
        }
        
        String[][] info = new String[4][pathwayTitles.size()];
        info[0] = pathwayTitles.toArray(info[0]);
        info[1] = pathwayHtmlFileNames.toArray(info[1]);
        info[2] = pathwayGpmlFileNames.toArray(info[2]);
        info[3] = pathwayPValues.toArray(info[3]);
        return info;
    }
    
    /**
     * Check if this pathway has no genes from the sample associated with it.
     * @return true if there are genes, false otherwise
     */
    public boolean noGenes() {
        if (this.genes.size()==0) {
            return true;
        }
        return false;
    }
    
    private String genesToString() {
        Iterator it = genes.iterator();
        StringBuffer genesString = new StringBuffer();
        if (it.hasNext()) {
            genesString.append(it.next());
        }
        while (it.hasNext()) {
            genesString.append(", "+it.next());
        }
        return genesString.toString();
    }
    
    public String getTitle() {
        return this.name;
    }
    
    public String getHtml() {
        return this.htmlFilePath;
    }
    
    public String getGpml() {
        return this.gpmlFilePath;
    }
    
    public String getPValue() {
        return String.format("%1.4e",this.p);
    }
    
    public double getQ() {
        return this.q;
    }
    
    private void setQ(double q) {
        this.q = q;
    }
    
    public double getCorrectedPValue(int numTests) {
        double correctedP = this.p * ((double) numTests);
        return correctedP;
    }
    
    public String toString() {
        return this.name;
    }
    
    public static String[] toStringArray(List<TestedPathway> testedPathways) {
        Iterator it = testedPathways.iterator();
        String[] pathwayArray = new String[testedPathways.size()];
        int counter = 0;
        while (it.hasNext()) {
            pathwayArray[counter++] = it.next().toString();
        }
        return pathwayArray;
    }
}

