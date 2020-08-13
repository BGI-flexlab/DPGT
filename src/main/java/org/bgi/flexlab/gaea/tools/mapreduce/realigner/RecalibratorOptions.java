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
package org.bgi.flexlab.gaea.tools.mapreduce.realigner;

import org.apache.commons.cli.ParseException;
import org.apache.hadoop.conf.Configuration;
import org.bgi.flexlab.gaea.data.mapreduce.options.HadoopOptions;
import org.bgi.flexlab.gaea.data.options.GaeaOptions;
import org.bgi.flexlab.gaea.tools.recalibrator.RecalibratorUtil;
import org.bgi.flexlab.gaea.tools.recalibrator.RecalibratorUtil.SolidNocallStrategy;
import org.bgi.flexlab.gaea.tools.recalibrator.RecalibratorUtil.SolidRecallMode;
import org.bgi.flexlab.gaea.tools.recalibrator.report.RecalibratorReportTable;
import org.bgi.flexlab.gaea.util.QualityUtils;

import java.util.ArrayList;
import java.util.List;

public class RecalibratorOptions extends GaeaOptions implements HadoopOptions {
	private final static String SOFTWARE_NAME = "BaseRecalibration";
	private final static String SOFTWARE_VERSION = "1.0";

	public RecalibratorOptions() {
		addOption("a", "defaultPlatform", true, "If a read has no platform then default to the provided String."
				+ " Valid options are illumina, 454, and solid.");
		addOption("A", "covariates", true,
				"One or more(separated by comma) covariates to be used in the recalibration.");
		addOption("B", "knownSites", true,
				"know  variant site,for example dbsnp!more than one file ,separation by Comma(,)!!");
		addOption("b", "bintag", true, "the binary tag covariate name if using it");
		addOption("C", "CachedRef", false, "cache reference");
		addOption("D", "ddq", true, "default quality for the base deletions covariate.(default:45)");
		addOption("E", "mdq", true, "default quality for the base mismatches covariate.(default:-1)");
		addOption("F", "idq", true, "default quality for the base insertions covariate.（default:45)");
		addOption("f", "forcePlatform", true, "If provided, the platform of EVERY read will be forced to be the "
				+ "provided String. Valid options are illumina, 454, and solid.");
		addOption("g", "mcs", true, "size of the k-mer context to be used for base mismatches.(default:2)");
		addOption("G", "ics", true,
				"size of the k-mer context to be used for base insertions and deletions.(default:3)");
		addOption("j", "recalibrationReference", true, "bqsr reference");
		addOption("k", "knownSite", true, "know  variant site,for example dbsnp!");
		addOption("K", "ls", false, "If specified, just list the available covariates and exit.");
		addOption("N", "noStandard", false,
				"If specified, do not use the standard set of covariates, but rather just the "
						+ "ones listed using the -cov argument.");
		addOption("p", "muq", true, "minimum quality for the bases to be preserved.");
		addOption("P", "sMode", true,
				"How should we recalibrate solid bases in which the reference was inserted? Options = DO_NOTHING, SET_Q_ZERO, "
						+ "SET_Q_ZERO_BASE_N, or REMOVE_REF_BIAS.");
		addOption("Q", "ql", true, "number of distinct quality scores in the quantized output.(default:16)");
		addOption("r", "reference", true, "bqsr reference");
		addOption("S", "solid_nocall_strategy", true,
				"Defines the behavior of the recalibrator when it encounters no calls in the color space. "
						+ "Options = THROW_EXCEPTION, LEAVE_READ_UNRECALIBRATED, or PURGE_READ.");
		addOption("T", "lqt", true,
				"minimum quality for the bases in the tail of the reads to be considered.(default:2)");
		addOption("w", "keyWindow", true, "window size for key[10000]");
		FormatHelpInfo(SOFTWARE_NAME, SOFTWARE_VERSION);
	}

	public boolean LIST_ONLY = false;

	public String[] COVARIATES = null;

	public boolean DO_NOT_USE_STANDARD_COVARIATES = false;

	public SolidRecallMode SOLID_RECAL_MODE = SolidRecallMode.SET_Q_ZERO;

	public SolidNocallStrategy SOLID_NOCALL_STRATEGY = SolidNocallStrategy.THROW_EXCEPTION;

	public int MISMATCHES_CONTEXT_SIZE;

	public int INDELS_CONTEXT_SIZE;

	public byte MISMATCHES_DEFAULT_QUALITY;

	public byte INSERTIONS_DEFAULT_QUALITY;

	public byte DELETIONS_DEFAULT_QUALITY;

	public byte LOW_QUALITY_TAIL;

	public int QUANTIZING_LEVELS;

	public int PRESERVE_QSCORES_LESS_THAN;

	public String BINARY_TAG_NAME = null;

	public String DEFAULT_PLATFORM = null;

	public String FORCE_PLATFORM = null;

	public boolean KEEP_INTERMEDIATE_FILES = false;

	public boolean NO_PLOTS = false;

	private List<String> knownSites = null;

	private String reference;

	private boolean isCachedRef;

	private int winSize;

	@Override
	public void setHadoopConf(String[] args, Configuration conf) {
		conf.setStrings("args", args);
	}

	@Override
	public void getOptionsFromHadoopConf(Configuration conf) {
		String[] args = conf.getStrings("args");
		this.parse(args);
	}

	public void parse(RecalibratorReportTable argsTable) {
		for (int i = 0; i < argsTable.getRowNumber(); i++) {
			final String argument = argsTable.get(i, "Argument").toString();
			Object value = argsTable.get(i, RecalibratorUtil.ARGUMENT_VALUE_COLUMN_NAME);
			if (value.equals("null"))
				value = null;

			if (argument.equals("covariate") && value != null)
				COVARIATES = value.toString().split(",");
			else if (argument.equals("no_standard_covs"))
				DO_NOT_USE_STANDARD_COVARIATES = Boolean.parseBoolean((String) value);
			else if (argument.equals("solid_recal_mode"))
				SOLID_RECAL_MODE = RecalibratorUtil.SolidRecallMode.recalModeFromString((String) value);
			else if (argument.equals("solid_nocall_strategy"))
				SOLID_NOCALL_STRATEGY = RecalibratorUtil.SolidNocallStrategy.nocallStrategyFromString((String) value);
			else if (argument.equals("mismatches_context_size"))
				MISMATCHES_CONTEXT_SIZE = Integer.parseInt((String) value);
			else if (argument.equals("indels_context_size"))
				INDELS_CONTEXT_SIZE = Integer.parseInt((String) value);
			else if (argument.equals("mismatches_default_quality"))
				MISMATCHES_DEFAULT_QUALITY = Byte.parseByte((String) value);
			else if (argument.equals("insertions_default_quality"))
				INSERTIONS_DEFAULT_QUALITY = Byte.parseByte((String) value);
			else if (argument.equals("deletions_default_quality"))
				DELETIONS_DEFAULT_QUALITY = Byte.parseByte((String) value);
			else if (argument.equals("low_quality_tail"))
				LOW_QUALITY_TAIL = Byte.parseByte((String) value);
			else if (argument.equals("default_platform"))
				DEFAULT_PLATFORM = (String) value;
			else if (argument.equals("force_platform"))
				FORCE_PLATFORM = (String) value;
			else if (argument.equals("quantizing_levels"))
				QUANTIZING_LEVELS = Integer.parseInt((String) value);
			else if (argument.equals("keep_intermediate_files"))
				KEEP_INTERMEDIATE_FILES = Boolean.parseBoolean((String) value);
			else if (argument.equals("no_plots"))
				NO_PLOTS = Boolean.parseBoolean((String) value);
			else if (argument.equals("binary_tag_name"))
				BINARY_TAG_NAME = (value == null) ? null : (String) value;
		}
	}

	@Override
	public void parse(String[] args) {
		try {
			cmdLine = parser.parse(options, args);
		} catch (ParseException e) {
			System.err.println(e.getMessage());
			printHelpInfotmation(SOFTWARE_NAME);
			System.exit(1);
		}

		setReference();
		setKnownSite();
		setSOLID_RECAL_MODE(getOptionValue("P", null));
		setSOLID_NOCALL_STRATEGY(getOptionValue("S", null));

		BINARY_TAG_NAME = getOptionValue("b", null);
		DEFAULT_PLATFORM = getOptionValue("a", null);
		FORCE_PLATFORM = getOptionValue("f", null);
		COVARIATES = getOptionValue("A", null) == null ? null : getOptionValue("A", null).split(",");
		MISMATCHES_CONTEXT_SIZE = getOptionIntValue("g", 2);
		INDELS_CONTEXT_SIZE = getOptionIntValue("G", 3);

		MISMATCHES_DEFAULT_QUALITY = (byte) getOptionIntValue("E", -1);
		INSERTIONS_DEFAULT_QUALITY = (byte) getOptionIntValue("F", 45);
		DELETIONS_DEFAULT_QUALITY = (byte) getOptionIntValue("D", 45);
		LOW_QUALITY_TAIL = (byte) getOptionIntValue("T", 2);

		PRESERVE_QSCORES_LESS_THAN = getOptionIntValue("p", QualityUtils.MINIMUM_USABLE_QUALITY_SCORE);
		QUANTIZING_LEVELS = getOptionIntValue("Q", 16);
		LIST_ONLY = getOptionBooleanValue("K", false);
		DO_NOT_USE_STANDARD_COVARIATES = getOptionBooleanValue("N", false);

		isCachedRef = getOptionBooleanValue("C", false);
		winSize = getOptionIntValue("w", 10000);
	}

	private void setKnownSite() {
		if (knownSites == null)
			knownSites = new ArrayList<String>();
		if (getOptionValue("B", null) != null) {
			String[] str = getOptionValue("B", null).split(",");
			for (String s : str)
				knownSites.add(s);
		} else
			this.knownSites.add(getOptionValue("k", null));

		if (knownSites.size() < 1)
			throw new RuntimeException("must set known dbsnp!");
	}

	private void setReference() {
		reference = getOptionValue("j", null);
		if (reference == null)
			reference = getOptionValue("r", null);

		if (reference == null)
			throw new RuntimeException("must set reference path!");
	}

	public int getWindowsSize() {
		return this.winSize;
	}

	public String getReferenceSequencePath() {
		return this.reference;
	}

	public List<String> getKnowSite() {
		return knownSites;
	}

	public SolidRecallMode getSOLID_RECAL_MODE() {
		return SOLID_RECAL_MODE;
	}

	@SuppressWarnings("static-access")
	public void setSOLID_RECAL_MODE(String optionValue) {
		if (optionValue == null)
			SOLID_RECAL_MODE = SOLID_RECAL_MODE.SET_Q_ZERO;
		else
			SOLID_RECAL_MODE = SOLID_RECAL_MODE.recalModeFromString(optionValue);
	}

	public SolidNocallStrategy getSOLID_NOCALL_STRATEGY() {
		return SOLID_NOCALL_STRATEGY;
	}

	public void setSOLID_NOCALL_STRATEGY(String optionValue) {
		if (optionValue == null)
			SOLID_NOCALL_STRATEGY = SolidNocallStrategy.THROW_EXCEPTION;
		else
			SOLID_NOCALL_STRATEGY = SolidNocallStrategy.nocallStrategyFromString(optionValue);
	}

	public boolean isCachedRef() {
		return isCachedRef;
	}
}
