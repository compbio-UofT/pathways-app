pathways-app
============

Pathways app for MedSavant.

##Install instructions

This app is available in the MedSavant app store.

If you would like to use the MedSavant Pathways plugin available in this repository, you need to have the MedSavant Client and Server already installed and functional. See http://genomesavant.com/p/medsavant/ for details.

Once you have the Client and Server ready, download MedSavant-App - Pathways-1.0.0.jar from the dist folder of this repository, and put in in your ~/.medsavant/plugins/ directory.

The next time you log into the MedSavant Client, the icon for Pathway Analysis should be displayed. Clicking on that icon will allow you to use the Pathways plugin.

##Usage instructions

A tutorial for the pathways plugin can be found [here](getting_started_guide.pdf).

##Build instructions

To build the project from source, you must have NetBeans installed. Download the source code, and go into the NetBeans menu to File > Open Project, and then navigate to the downloaded code. Select the project in the Project panel, and select Run > Clean and Build.

##Future Directions
Features that would be nice to have in the future include:
* Ability to print out a comprehensive PDF report of results
* [GREAT annotation](http://bejerano.stanford.edu/papers/GREAT.pdf) for intergenic variants
* BiNGO test for Gene Ontology term enrichment
* incorporate functional impact scores into enrichment tests
