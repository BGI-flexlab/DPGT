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
package org.bgi.flexlab.gaea.data.options;

import org.apache.commons.cli.*;

import java.util.ArrayList;
import java.util.List;

public abstract class GaeaOptions {
	protected Options options = new Options();
	protected CommandLine cmdLine;
	protected CommandLineParser parser = new PosixParser();
	protected HelpFormatter helpInfo = new HelpFormatter();
	protected List<Option> optionsList = new ArrayList<Option>();
	
	/**
	 * parse parameter
	 */
	abstract public void parse(String[] args);

	/**
	 * help info must use it before parse.
	 */
	public void FormatHelpInfo(String softwareName, String version) {
		StringBuilder sb = new StringBuilder();
		if (softwareName != null) {
			sb.append("Software name: ");
			sb.append(softwareName);
		}
		sb.append("\nVersion: ");
		sb.append(version);
		sb.append("\nLast update: 2020.08.11\n");
		sb.append("Developed by: BGI-shenzhen\n");
		sb.append("Authors: Gong Chun\n");
		sb.append("E-mail: gongchun@genomics.org.cn\n");
		sb.append("Copyright(c) 2020: BGI. All Rights Reserved.\n\n");
		helpInfo.setNewLine("\n");
		if(softwareName == null)
			softwareName = "tools_name";
		helpInfo.setSyntaxPrefix("hadoop jar " + "DGPT.jar " + softwareName
				+ " [options]\n" + sb.toString() + "\n");
		helpInfo.setWidth(2 * HelpFormatter.DEFAULT_WIDTH);
	}
	
	protected void initialization() {
		
	}

	protected void printHelpInfotmation(String softwareName) {
		helpInfo.printHelp(softwareName + " options list:", options);
	}

	/**
	 * add boolean option
	 */
	protected void addBooleanOption(String opt, String longOpt,
			String description) {
		addOption(opt, longOpt, false, description);
	}

	/**
	 * add normal option.
	 * 
	 * @param opt
	 * @param longOpt
	 * @param hasArg
	 * @param description
	 */
	protected void addOption(String opt, String longOpt, boolean hasArg,
			String description) {
		addOption(opt, longOpt, hasArg, description, false);
	}

	/**
	 * add required option
	 */
	protected void addOption(String opt, String longOpt, boolean hasArg,
			String description, boolean required) {
		Option option = new Option(opt, longOpt, hasArg, description);
		option.setRequired(required);
		addOption(option);
	}
	
	protected void addOption(Option option){
		options.addOption(option);
		optionsList.add(option);
	}

	protected String getOptionValue(String opt, String defaultValue) {
		if (cmdLine.hasOption(opt))
			return cmdLine.getOptionValue(opt);

		return defaultValue;
	}

	protected int getOptionIntValue(String opt, int defaultValue) {
		if (cmdLine.hasOption(opt))
			return Integer.parseInt(cmdLine.getOptionValue(opt));
		return defaultValue;
	}

	protected boolean getOptionBooleanValue(String opt, boolean defaultValue) {
		if (cmdLine.hasOption(opt))
			return true;
		return defaultValue;
	}

	protected double getOptionDoubleValue(String opt, double defaultValue) {
		if (cmdLine.hasOption(opt))
			return Double.parseDouble(cmdLine.getOptionValue(opt));
		return defaultValue;
	}

	protected long getOptionLongValue(String opt, long defaultValue) {
		if (cmdLine.hasOption(opt))
			return Long.parseLong(cmdLine.getOptionValue(opt));
		return defaultValue;
	}

	protected byte getOptionByteValue(String opt, byte defaultValue) {
		if (cmdLine.hasOption(opt))
			return Byte.parseByte(cmdLine.getOptionValue(opt));
		return defaultValue;
	}

	protected short getOptionShortValue(String opt, short defaultValue) {
		if (cmdLine.hasOption(opt))
			return Short.parseShort(cmdLine.getOptionValue(opt));
		return defaultValue;
	}
	
	public Options getOptions(){
		return this.options;
	}
	
	public List<Option> getOptionList(){
		return this.optionsList;
	}
}
