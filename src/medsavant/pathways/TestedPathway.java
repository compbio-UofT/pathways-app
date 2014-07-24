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
    private double p;
    private String name;
    private HashSet<String> genes;
    private String htmlFilePath;
    private String pngFilePath;
    private String gpmlFilePath;
    private int totalGenes;
    public static int PATHWAYNAMEINDEX = 1;
    public static int PVALUEINDEX = 0;
    public static int GENESINDEX = 2;
    public static int NUMGENESINDEX = 3;
    private double q;
    /**
     * Initialize tested pathway.
     * @param p double p-value of tested pathway from hypergeometric test results.
     * @param name pathway title
     * @param genes set of genes involved in pathway
     */
    /*public TestedPathway(double p, String name, HashSet<String> genes) {
        this.p = p;
        this.name = name;
        this.genes = new HashSet<String>(genes);
    }*/
    
    
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
    
    public Object[] getObjectList() {
        if (genes.size() == 0) {
            return null;
        }
        Object[] objectArray = new Object[4];
        objectArray[0] = this.getPValue();
        objectArray[1] = this.name;
        objectArray[2] = genesToString();
        objectArray[3] = this.totalGenes;
        return objectArray;
    }
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
    
    private static void assignQValues(List<TestedPathway> list, double fdrCutoff) {
        Collections.sort(list, new PathwayEnrichmentProbabilityComparator());
        int size = list.size();
        for (int i = 0; i < size; i++) {
            list.get(i).setQ(((double) i)*fdrCutoff/((double) size));
        }
    }
    private Object[] getBHObjectList() {
        if (genes.size() == 0) {
            return null;
        }
        Object[] objectArray = new Object[4];
        objectArray[0] = this.getQ();
        objectArray[1] = this.name;
        objectArray[2] = genesToString();
        objectArray[3] = this.totalGenes;
        return objectArray;
    }
    public static List<Object[]> convertToObjectListBH(List<TestedPathway> list, double fdrCutoff) {
        // MAKE BENJAMINI HOCHBERG P-VALUES
        assignQValues(list, fdrCutoff);
        List<Object[]> newList = new ArrayList<Object[]>();
        Iterator it = list.iterator();
        Object[] objectArray;
        int numTests = list.size();
        while (it.hasNext()) {
            objectArray = ((TestedPathway) it.next()).getBHObjectList();
            if (objectArray!=null) {
                newList.add(objectArray);
                
            }
        }
        return newList;
    }
    
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
                System.out.println("pathway title: "+pathway.getTitle()+" pathway html: "+pathway.getHtml()+" pathway gpml:"+pathway.getGpml());
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
    public String getCorrectedPValue(int numTests) {
        double correctedP = this.p * ((double) numTests);
        return String.format("%1.4e",correctedP);
    }
}

