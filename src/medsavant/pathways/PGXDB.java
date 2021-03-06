package medsavant.pathways;


import java.io.File;
import java.io.FileInputStream;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.ResultSetMetaData;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.LinkedList;
import java.util.List;
import org.ut.biolab.medsavant.client.settings.DirectorySettings;


/**
 * Local pharmacogenomics DB.
 * 
 * @author rammar
 */
public class PGXDB {
	
	private static final String db_path_prefix= "mem:pharmacogenomicsdb";
	private static final String DB_URL= "jdbc:hsqldb:" + db_path_prefix;
	private static final String GENE_MARKER_LIST_FILE_PATH= "/pgx/localDBFiles/gene_marker_list.txt";
	private static final String HAPLOTYPE_MARKERS_FILE_PATH= "/pgx/localDBFiles/haplotype_markers.txt";
	private static final String MARKER_COORDINATES_FILE_PATH= "/pgx/localDBFiles/marker_coordinates.txt";
	private static final String HAPLOTYPE_ACTIVITY_FILE_PATH= "/pgx/localDBFiles/haplotype_activity.txt";
	private static final String ACTIVITY_TO_METABOLIZER_FILE_PATH= "/pgx/localDBFiles/activity_to_metabolizer.txt";
	private static final String PHENOTYPE_TO_METABOLIZER_FILE_PATH= "/pgx/localDBFiles/phenotype_to_metabolizer.txt";
	private static final String GENE_MARKER_LIST_TABLE_NAME= "gene_marker_list";
	private static final String HAPLOTYPE_MARKERS_TABLE_NAME= "haplotype_markers";
	private static final String MARKER_COORDINATES_TABLE_NAME= "marker_coordinates";
	private static final String HAPLOTYPE_ACTIVITY_TABLE_NAME= "haplotype_activity";
	private static final String ACTIVITY_TO_METABOLIZER_TABLE_NAME= "activity_to_metabolizer";
	private static final String PHENOTYPE_TO_METABOLIZER_TABLE_NAME= "phenotype_to_metabolizer";
	
	private static Connection conn;
	
			
	/**
	 * Creates all tables to the DB and populates with all data from packaged
	 * text files.
	 * @throws SQLException 
	 */
	public static void initialize() throws SQLException {		
		conn= connectionToServer();	
		createSchema(conn);
		loadTables(conn);
	}

	
	/** 
	 * Close connection to DB.
	 */
	public static void closeConnectionToDB() throws SQLException {
		conn.close();
	}
	
	
	/**
	 * Create JDBC-based connection to HSQL server.
	 * By default, connection is auto-committing.
	 * @return the Connection to the local HSQLDB JDBC driver
	 */
	private static Connection connectionToServer() throws SQLException {
		/* Snippet of code from HyperSQL documentation. Loads the HSQLDB JDBC driver. */
		try {
			Class.forName("org.hsqldb.jdbcDriver");
		} catch (Exception e) {
			System.err.println("[" + PGXDB.class.getSimpleName() + 
				"]: failed to load HSQLDB JDBC driver.");
			e.printStackTrace();
		}
		
		return DriverManager.getConnection(DB_URL);
	}
	
	
	/**
	 * Creates the table schema.
	 * @param c The HSQL DB connection
	 */
	private static void createSchema(Connection c) {
		try {
			Statement s= c.createStatement();
			String sql;

			/* Execute SQL as a batch, rather than individually. This is more efficient
			 * when processing many. */
			
			
			/* Enable case insensitivity. According to the hsqldb documentation,
			 * this must be switched before creating tables:
			 * http://www.hsqldb.org/doc/guide/ch09.html#set_ignorecase-section */
			sql=	"SET IGNORECASE TRUE";
			s.addBatch(sql);
			
			/* Create all the tables. These will be populated by CSV Loader. */
			
			sql=	"CREATE TABLE " + GENE_MARKER_LIST_TABLE_NAME + " ( " +
					"  Gene varchar(20) NOT NULL, " +
					"  Marker_list varchar(50000) NOT NULL, " + // could be a very long field
					"  PRIMARY KEY (Gene)  " +
					")";
			s.addBatch(sql);

			
			/* By default *1 (no markers in the marker_info column) is an empty
			 * string, rather than NULL */
			sql=	"CREATE TABLE " + HAPLOTYPE_MARKERS_TABLE_NAME + " ( " +
					"	Gene varchar(20) NOT NULL, " + 
					"	Haplotype_ID varchar(1000) DEFAULT NULL, " +
					"	Haplotype_Symbol varchar(100) NOT NULL, " + 
					"	Marker_info varchar(50000) DEFAULT NULL, " + // could be a very long field
					"	PRIMARY KEY (Gene,Haplotype_Symbol)  " +
					")";
			s.addBatch(sql);
			
			sql=	"CREATE TABLE " + MARKER_COORDINATES_TABLE_NAME + " ( " +
					"	Marker varchar(40) NOT NULL, " + 
					"	Chromosome varchar(20) NOT NULL, " +
					"	Position int NOT NULL, " +
					"	Ref varchar(10000) NOT NULL, " + // ref and alt can be very long nucleotide strings
					"	Alt varchar(10000) NOT NULL, " +			
					"	PRIMARY KEY (Marker,Position, Alt)  " + // same rsID can indicate multiple alt calls and have multiple positions (see rs72558184)
					")";
			s.addBatch(sql);
			
			sql=	"CREATE TABLE " + HAPLOTYPE_ACTIVITY_TABLE_NAME + " ( " +
					"	Gene varchar(20) NOT NULL, " +
					"	Haplotype varchar(10) NOT NULL, " +
					"	Activity_Score double DEFAULT NULL, " +
					"	Activity_Phenotype varchar(100) NOT NULL, " +
					"	Pubmed_ID varchar(1000) NOT NULL, " +
					"	PRIMARY KEY (Gene, Haplotype) " +
					")";
			s.addBatch(sql);

			sql=	"CREATE TABLE " + ACTIVITY_TO_METABOLIZER_TABLE_NAME + " ( " +
					"	Total_Activity_Score_Minimum double NOT NULL, " +
					"	Metabolizer_class varchar(100) NOT NULL, " +
					")";
			s.addBatch(sql);
			
			sql=	"CREATE TABLE " + PHENOTYPE_TO_METABOLIZER_TABLE_NAME + " ( " +
					"	Haplotype_1_activity varchar(100) NOT NULL, " +
					"	Haplotype_2_activity varchar(100) NOT NULL, " +
					"	Metabolizer_class varchar(100) NOT NULL, " +
					")";
			s.addBatch(sql);
			
			s.executeBatch();
			s.close();
		} catch (SQLException e) {
			System.err.println("[" + PGXDB.class.getSimpleName() + "]: Error creating tables " + e.toString());
			e.printStackTrace();
		}
	}
	
	
	/** 
	 * Load tables into the DB.
	 * @param c The HSQL DB connection
	 */
	private static void loadTables(Connection c) throws SQLException {		
		/* Load the delimited tables from text files. */
		CSVLoader loader;
		try {			
			loader= new CSVLoader(connectionToServer()); // pass a new connection since it auto-closes it.
			loader.setSeprator('\t');
			loader.loadCSV(PGXDB.class.getResourceAsStream(GENE_MARKER_LIST_FILE_PATH),
				GENE_MARKER_LIST_TABLE_NAME, false);
			
			loader= new CSVLoader(connectionToServer()); // pass a new connection since it auto-closes it.
			loader.setSeprator('\t');
			loader.loadCSV(PGXDB.class.getResourceAsStream(HAPLOTYPE_MARKERS_FILE_PATH),
				HAPLOTYPE_MARKERS_TABLE_NAME, false);
			
			loader= new CSVLoader(connectionToServer()); // pass a new connection since it auto-closes it.
			loader.setSeprator('\t');
			loader.loadCSV(PGXDB.class.getResourceAsStream(MARKER_COORDINATES_FILE_PATH),
				MARKER_COORDINATES_TABLE_NAME, false);
			
			loader= new CSVLoader(connectionToServer()); // pass a new connection since it auto-closes it.
			loader.setSeprator('\t');
			loader.loadCSV(PGXDB.class.getResourceAsStream(HAPLOTYPE_ACTIVITY_FILE_PATH),
				HAPLOTYPE_ACTIVITY_TABLE_NAME, false);
			
			loader= new CSVLoader(connectionToServer()); // pass a new connection since it auto-closes it.
			loader.setSeprator('\t');
			loader.loadCSV(PGXDB.class.getResourceAsStream(ACTIVITY_TO_METABOLIZER_FILE_PATH),
				ACTIVITY_TO_METABOLIZER_TABLE_NAME, false);
			
			loader= new CSVLoader(connectionToServer()); // pass a new connection since it auto-closes it.
			loader.setSeprator('\t');
			loader.loadCSV(PGXDB.class.getResourceAsStream(PHENOTYPE_TO_METABOLIZER_FILE_PATH),
				PHENOTYPE_TO_METABOLIZER_TABLE_NAME, false);
			
		} catch (Exception e) {
			System.err.println("[" + PGXDB.class.getSimpleName() + "]: Error loading tables " + e.toString());
			e.printStackTrace();
		}
	}
	
	
	/**
	 * Return this database's Connection.
	 * @return this database's Connection; null if connection has not been initialized.
	 */
	public static Connection getConnection() {
		return conn;
	}
	
	
	/**
	 * Execute an SQL command.
	 * @param	sql	The SQL statement to execute.
	 * @return	ResultSet object corresponding to the SQL statement's output
	 * @precondition Static Connection conn not null.
	 */
	public static ResultSet executeQuery(String sql) throws SQLException {
		Statement s= conn.createStatement();
		s.execute(sql);
		
		ResultSet rs= s.getResultSet();
		
		s.close();
		
		return rs;
	}
	
	
	/**
	 * Print text output for a given query.
	 * @param	sql	String SQL statement to be executed. Includes all possible SQL statements.
	 * @precondition Static Connection conn not null.
	 */
	public static void printQueryResults(String sql) throws SQLException {
		Statement s= conn.createStatement();
		s.execute(sql);
		
		ResultSet rs= s.getResultSet();
		while (rs.next()) {
			System.out.println("[" + PGXDB.class.getSimpleName() + "]: " + getRowAsString(rs));
		}
		
		s.close();
	}
	
	
	/**
	 * Output the entire row because ResultSet get methods are busting my balls.
	 */
	public static String getRowAsString(ResultSet rs) throws SQLException {
		String result= "";
		
		ResultSetMetaData rsmd= rs.getMetaData();
		int columnsNumber= rsmd.getColumnCount();
		
		// Remember, columsn start at index 1
		for (int i= 1; i <= columnsNumber; ++i) {
			result += rs.getString(i) + "\t";
		}
		
		return result;
	}
	
	
	/**
	 * Output the entire row because ResultSet get methods are busting my balls.
	 * Gets the row pointed to by cursor's current position.
	 */
	public static List<Object> getRowAsList(ResultSet rs) throws SQLException {
		ResultSetMetaData rsmd= rs.getMetaData();
		int columnsNumber= rsmd.getColumnCount();
		
		LinkedList<Object> result= new LinkedList<Object>();
		
		// Remember, columns start at index 1
		for (int i= 1; i <= columnsNumber; ++i) {
			result.add(rs.getString(i));
		}
		
		return result;
	}
	
}
