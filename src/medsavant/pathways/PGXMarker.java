/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package medsavant.pathways;

	/**
	 * Inner class to represent PGx markers.
	 */
	public class PGXMarker {
		public String markerID;
		public String chromosome;
		public String position;
		public String ref;
		public String alt;
		
		public PGXMarker(String markerID, String chromosome, String position, String ref, String alt) {
			this.markerID= markerID;
			this.chromosome= chromosome;
			this.position= position;
			this.ref= ref;
			this.alt= alt;
		}
	}