/*******************************************************************************
 * Copyright (c) 2017, BGI-Shenzhen
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 *******************************************************************************/
package org.bgi.flexlab.gaea.tools.annotator.config;

import com.google.gson.Gson;
import com.google.gson.JsonSyntaxException;
import com.google.gson.reflect.TypeToken;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.bgi.flexlab.gaea.data.structure.reference.ReferenceShare;
import org.bgi.flexlab.gaea.tools.annotator.codons.CodonTable;
import org.bgi.flexlab.gaea.tools.annotator.codons.CodonTables;
import org.bgi.flexlab.gaea.tools.annotator.effect.SnpEffectPredictor;
import org.bgi.flexlab.gaea.tools.annotator.interval.Chromosome;
import org.bgi.flexlab.gaea.tools.annotator.interval.Genome;
import org.bgi.flexlab.gaea.tools.annotator.util.CountByType;
import org.bgi.flexlab.gaea.tools.annotator.util.Gpr;
import org.bgi.flexlab.gaea.tools.annotator.util.Timer;
import org.bgi.flexlab.gaea.tools.mapreduce.annotator.AnnotatorOptions;
import org.bgi.flexlab.gaea.util.OrderedProperties;

import java.io.*;
import java.util.*;

public class Config implements Serializable {
	
	private static final long serialVersionUID = 2793968455002586646L;
	
	private static Config configInstance = null; 
	
	public static final String KEY_REFERENCE = "ref";
	public static final String KEY_GENEINFO_PREFIX = "GeneInfo";
	public static final String KEY_VARIANT_PREFIX = "Variant";
	public static final String KEY_EFFECT_NOMEN = "EffectNomen";
	public static final String KEY_DEFAULT_VALUE = "SetDefaultValue";
	public static final String KEY_TSV_PREFIX = "TSV";
	public static final String KEY_CODON_PREFIX = "codon.";
	public static final String KEY_CODONTABLE_SUFIX = ".codonTable";
	public static final String DB_CONFIG_JSON = "AnnotatorConfig.json";
	public static final String ANNO_FIELDS_SUFIX = ".fields";
	public static int MAX_WARNING_COUNT = 20;
	
	private String  ref = null;
	private String  geneInfo = null;

	private boolean debug = false; // Debug mode?
	private boolean verbose = false; // Verbose
	private boolean treatAllAsProteinCoding;
	private boolean onlyRegulation; // Only use regulation features
	private boolean errorOnMissingChromo; // Error if chromosome is missing
	private boolean errorChromoHit; // Error if chromosome is not hit in a query
	private boolean hgvs = true; // Use HGVS notation?
	private boolean hgvsShift = true; // Shift variants according to HGVS notation (towards the most 3prime possible coordinate)
	private boolean hgvsOneLetterAa = false; // Use HGVS 1 letter amino acid in HGVS notation?
	private boolean hgvsTrId = false; // Use HGVS transcript ID in HGVS notation?
	private CountByType warningsCounter = new CountByType();
	private HashMap<String, ArrayList<String>> annoFieldsByDB = null;  // fields in user config
	private DatabaseJson databaseJson;
	private List<String> dbNameList;
	private List<String> fields = new ArrayList<>();
	private List<String> fieldsWithoutVariant = null;
	private Map<String, String> headerByField = new HashMap<>();
	private Map<String, String> defaultValue = new HashMap<>();
	private Properties properties;
	private Genome genome;
	private SnpEffectPredictor snpEffectPredictor;
	private Configuration conf;
	private AnnotatorOptions options;
    private boolean useSimpleEffectNom;

    public Config(Configuration conf) throws IOException {
		this.conf = conf;
		init();
		configInstance = this;
	}
	
	public Config(Configuration conf, ReferenceShare genomeShare) throws IOException {
		this.conf = conf;
		init();
		configInstance = this;
		genome = new Genome(ref,genomeShare);
	}
	
	private void init() throws IOException {
		options = new AnnotatorOptions();
		options.getOptionsFromHadoopConf(conf);
		treatAllAsProteinCoding = false;
		onlyRegulation = false;
		useSimpleEffectNom = false;
		errorOnMissingChromo = true;
		errorChromoHit = true;
		verbose = options.isVerbose();
		debug = options.isDebug();

		loadProperties(options.getConfigFile()); // Read config file and get a genome
//		TODO 支持在配置文件中自定义密码子体系 - CodonTable
//		createCodonTables(genomeVersion, properties);  
		parseProperties();
	}

	/**
	 * Load properties from configuration file
	 * @return true if success
	 */
	boolean loadProperties(String configFileName) {
		properties = new OrderedProperties();
		try {
			Path confFilePath = new Path(configFileName);
			FileSystem fs = confFilePath.getFileSystem(conf);
			if(!fs.exists(confFilePath)) {
				throw new RuntimeException(confFilePath.toString() + " don't exist.");
			}
			if(!fs.isFile(confFilePath)) {
				throw new RuntimeException(confFilePath.toString() + " is not a file.");
			}
			properties.load(fs.open(confFilePath));
			
			if (!properties.isEmpty()) {
				return loadJson();
			}
		} catch (Exception e) {
			properties = null;
			throw new RuntimeException(e);
		}

		return false;
	}
	
	/**
	 * load database info json
	 * @return
	 */
	private boolean loadJson() {

		Gson gson = new Gson();
		try {
			Reader reader;
			if(properties.containsKey("DB_CONFIG")){
				String databaseConfig = properties.getProperty("DB_CONFIG");
				reader = new InputStreamReader(new FileInputStream(databaseConfig));
			}else {
				reader = new InputStreamReader(Config.class.getClassLoader().getResourceAsStream(DB_CONFIG_JSON), "UTF-8");
			}
			databaseJson = gson.fromJson(reader, new TypeToken<DatabaseJson>(){}.getType());
		} catch ( JsonSyntaxException | IOException e) {
			databaseJson = null;
			Timer.showStdErr("Read a wrong json config file:" + DB_CONFIG_JSON);
			e.printStackTrace();
		}
		return databaseJson != null;
	}
	
	/**
	 * Extract and create codon tables
	 */
	void createCodonTables(String genomeId, Properties properties) {
		//---
		// Read codon tables
		//---
		for (Object key : properties.keySet()) {
			if (key.toString().startsWith(KEY_CODON_PREFIX)) {
				String name = key.toString().substring(KEY_CODON_PREFIX.length());
				String table = properties.getProperty(key.toString());
				CodonTable codonTable = new CodonTable(name, table);
				CodonTables.getInstance().add(codonTable);
			}
		}

		//---
		// Assign codon tables for different genome+chromosome
		//---
		for (Object key : properties.keySet()) {
			String keyStr = key.toString();
			if (keyStr.endsWith(KEY_CODONTABLE_SUFIX) && keyStr.startsWith(genomeId + ".")) {
				// Everything between gneomeName and ".codonTable" is assumed to be chromosome name
				int chrNameEnd = keyStr.length() - KEY_CODONTABLE_SUFIX.length();
				int chrNameStart = genomeId.length() + 1;
				int chrNameLen = chrNameEnd - chrNameStart;
				String chromo = null;
				if (chrNameLen > 0) chromo = keyStr.substring(chrNameStart, chrNameEnd);

				// Find codon table
				String codonTableName = properties.getProperty(key.toString());
				CodonTable codonTable = CodonTables.getInstance().getTable(codonTableName);
				if (codonTable == null) throw new RuntimeException("Error parsing property '" + key + "'. No such codon table '" + codonTableName + "'");

				if (chromo != null) {
					// Find chromosome
					Chromosome chr = genome.getOrCreateChromosome(chromo);
					CodonTables.getInstance().set(genome, chr, codonTable);
				} else {
					// Set genome-wide chromosome table
					CodonTables.getInstance().set(genome, codonTable);
				}
			}
		}
	}
	
	private boolean parseProperties() throws IOException {
		
		// Sorted keys
		Set<String> keys = properties.stringPropertyNames();
		ref = properties.getProperty(KEY_REFERENCE);
		if(properties.containsKey(KEY_EFFECT_NOMEN))
			setUseSimpleEffectNom(properties.getProperty(KEY_EFFECT_NOMEN));
		setGeneInfo(properties.getProperty(KEY_GENEINFO_PREFIX));
		if(properties.containsKey(KEY_DEFAULT_VALUE))
			setDefaultValue(properties.getProperty(KEY_DEFAULT_VALUE));
		
		annoFieldsByDB = new HashMap<>();
		dbNameList = new ArrayList<>();

		//用户配置文件中注释字段的配置格式： dbName.N.fields = field1:header1,field2,field3
		for (String key : keys) {
			if(options.getOutputFormat() != AnnotatorOptions.OutputFormat.TSV && key.startsWith(KEY_TSV_PREFIX))
				continue;

			if (key.endsWith(ANNO_FIELDS_SUFIX)) {
				String dbName = key.substring(0, key.length() - ANNO_FIELDS_SUFIX.length());
				if(dbName.contains(".")){
					dbName = key.substring(0, dbName.lastIndexOf('.'));
				}
				String[] annoFields = properties.getProperty(key).split(",");
				ArrayList<String> annoFieldList = new ArrayList<>();

				for (String annoField : annoFields) {
					annoField = annoField.trim();
					if (annoField.contains(":")) {
						String[] tags = annoField.split(":");
						fields.add(tags[0]);
						headerByField.put(tags[0], tags[1]);
						annoFieldList.add(tags[0]);
					} else {
						annoFieldList.add(annoField);
						fields.add(annoField);
					}
				}

				if (!key.startsWith(KEY_GENEINFO_PREFIX) && !key.startsWith(KEY_VARIANT_PREFIX)) {
					if(!dbNameList.contains(dbName))
						dbNameList.add(dbName);
				}

				if(annoFieldsByDB.containsKey(dbName))
					annoFieldsByDB.get(dbName).addAll(annoFieldList);
				else
					annoFieldsByDB.put(dbName, annoFieldList);
			}
		}
		
		return true;
	}

	public void setUseSimpleEffectNom(String useSimpleEffectNomStr) {
		switch (useSimpleEffectNomStr.trim().toLowerCase()) {
			case "false":
				this.useSimpleEffectNom = false;
			case "f":
				this.useSimpleEffectNom = false;
			case "0":
				this.useSimpleEffectNom = false;
			default:
				this.useSimpleEffectNom = true;
		}
	}

	public Map<String, String> getDefaultValue() {
		return defaultValue;
	}

	public void setDefaultValue(String defaultValueStr) {
		String[] dv = defaultValueStr.split(",");
		for (String aDv : dv) {
			String[] kv = aDv.split(":", 2);
			if (kv.length == 2)
				defaultValue.put(kv[0], kv[1]);
		}
	}

	public boolean isUseSimpleEffectNom() {
		return useSimpleEffectNom;
	}


	public List<String> getFields() {
		return fields;
	}

	public List<String> getFieldsWithoutVariant() {
		if(fieldsWithoutVariant != null)
			return fieldsWithoutVariant;

		fieldsWithoutVariant = new ArrayList<>();
		for (String field: getFields())
			if(!annoFieldsByDB.get(KEY_VARIANT_PREFIX).contains(field))
				fieldsWithoutVariant.add(field);
		return fieldsWithoutVariant;
	}

	public Map<String, String> getHeaderByField() {
		return headerByField;
	}


	public String getHeaderNameByField(String field) {
		if(headerByField.containsKey(field))
			return headerByField.get(field);
		return field;
	}

	public ArrayList<String> getFieldsByDB(String dbName){
		return annoFieldsByDB.get(dbName);
	}

	public static Config getConfigInstance() {
		return configInstance;
	}

	public static Config get(){
		return configInstance;
	}

	public Genome getGenome() {
		return genome;
	}

	public String getGenomeVersion() {
		return ref;
	}
	
	public void setHgvsOneLetterAA(boolean hgvsOneLetterAa) {
		this.hgvsOneLetterAa = hgvsOneLetterAa;
	}

	public void setHgvsShift(boolean hgvsShift){
		this.hgvsShift = hgvsShift;
	}

	public void setHgvsTrId(boolean hgvsTrId) {
		this.hgvsTrId = hgvsTrId;
	}

	public void setOnlyRegulation(boolean onlyRegulation) {
		this.onlyRegulation = onlyRegulation;
	}

	public void setString(String propertyName, String value) {
		properties.setProperty(propertyName, value);
	}

	public void setTreatAllAsProteinCoding(boolean treatAllAsProteinCoding) {
		this.treatAllAsProteinCoding = treatAllAsProteinCoding;
	}

	public void setUseHgvs(boolean useHgvs) {
		hgvs = useHgvs;
	}

	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}
	
	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	public boolean isDebug() {
		return debug;
	}

	public boolean isErrorChromoHit() {
		return errorChromoHit;
	}

	public boolean isErrorOnMissingChromo() {
		return errorOnMissingChromo;
	}

	public boolean isHgvs() {
		return hgvs;
	}

	public boolean isHgvs1LetterAA() {
		return hgvsOneLetterAa;
	}

	public boolean isHgvsShift() {
		return hgvsShift;
	}

	public boolean isHgvsTrId() {
		return hgvsTrId;
	}

	public boolean isOnlyRegulation() {
		return onlyRegulation;
	}

	public boolean isTreatAllAsProteinCoding() {
		return treatAllAsProteinCoding;
	}

	public boolean isVerbose() {
		return verbose;
	}
	
	/**
	 * Show a warning message and exit
	 */
	public void warning(String warningType, String details) {
		long count = warningsCounter.inc(warningType);

		if (debug || count < MAX_WARNING_COUNT) {
			if (debug) Gpr.debug(warningType + details, 1);
			else System.err.println(warningType + details);
		} else if (count <= MAX_WARNING_COUNT) {
			String msg = "Too many '" + warningType + "' warnings, no further warnings will be shown:\n" + warningType + details;
			if (debug) Gpr.debug(msg, 1);
			else System.err.println(msg);
		}
	}

	public String getGeneInfo() {
		return geneInfo;
	}

	public void setGeneInfo(String geneInfo) {
		if (geneInfo.startsWith("/")) {
			this.geneInfo = "file://" + geneInfo;
		}else {
			this.geneInfo = geneInfo;
		}
	}

	public String getRef() {
		switch(ref){
		
			case "hg19":
			case "GRCh37":
				   return "GRCh37";
				   
			case "hg20":
			case "hg38":
			case "GRCh38":
				return "GRCh38";
			default:
					throw new RuntimeException("ref '" + ref + "' not found.");
		}
	}

	public void setRef(String ref) {
		this.ref = ref;
	}

	public SnpEffectPredictor getSnpEffectPredictor() {
		return snpEffectPredictor;
	}

	public void setSnpEffectPredictor(SnpEffectPredictor snpEffectPredictor) {
		this.snpEffectPredictor = snpEffectPredictor;
	}

	public List<String> getDbNameList() {
		return dbNameList;
	}

	public void setDbNameList(List<String> dbNameList) {
		this.dbNameList = dbNameList;
	}

	public DatabaseJson getDatabaseJson() {
		return databaseJson;
	}

	public String getHeaderString(){
		return getHeaderString("#", "\t");
	}

	public String getHeaderString(String prefix, String delimiter){
		ArrayList<String> headers = new ArrayList<>();
		for(String field: getFields()){
			headers.add(getHeaderNameByField(field));
		}
		return prefix+String.join(delimiter, headers);
	}

	public String getVCFHeaderString(){
		ArrayList<String> headers = new ArrayList<>();
		for(String field: getFieldsWithoutVariant()){
			headers.add(getHeaderNameByField(field));
		}
		return String.join("|", headers);
	}


	public AnnotatorOptions getOptions() {
		return options;
	}

	//	public static void main(String[] args) throws Exception {
//		Config config = new Config();
//		config.loadJson();
//		System.out.println(config.databaseJson.getDatabaseInfo("dbNSFP").getRefTable("GRCh38").getIndexTable());
//		System.out.println(config.databaseJson.getDatabaseInfo("dbNSFP").getFields());
//		if (config.databaseJson.getDatabaseInfo("dbNSFP").getFields() == null) {
//			System.out.println("dbNSFP fields is null!");
//		}else if (config.databaseJson.getDatabaseInfo("dbNSFP").getFields().isEmpty()) {
//			System.out.println("dbNSFP fields is empty!");
//		}
//		System.out.println("Good");
//	}

}
