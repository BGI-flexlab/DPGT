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
package org.bgi.flexlab.gaea.data.exception;

import htsjdk.samtools.SAMRecord;

import java.io.File;

public class UserException extends RuntimeException {
	private static final long serialVersionUID = 4451398435363715205L;

	public UserException(String msg) {
		super(msg);
	}

	public UserException(String msg, Throwable e) {
		super(msg, e);
	}

	protected static String getMessage(Throwable t) {
		String message = t.getMessage();
		return message != null ? message : t.getClass().getName();
	}
	
	public static class BadInput extends UserException {

		private static final long serialVersionUID = 1708634772791836896L;

		public BadInput(String message) {
            super(String.format("Bad input: %s", message));
        }
    }
	
	public static class PileupException extends UserException{

		private static final long serialVersionUID = 1L;

		public PileupException(String message) {
			super(String.format("pileup exception :", message));
		}
	}
	
	public static class MalformedVCFHeader extends UserException {
	        /**
			 * 
			 */
			private static final long serialVersionUID = 7027922936632516646L;

			public MalformedVCFHeader(String message) {
	            super(String.format("The provided VCF file has a malformed header: %s", message));
	        }
    }
	
	public static class MalformedVCF extends UserException {
	        /**
			 * 
			 */
			private static final long serialVersionUID = 8121781154975706475L;

			public MalformedVCF(String message, String line) {
	            super(String.format("The provided VCF file is malformed at line %s: %s", line, message));
	        }

	        public MalformedVCF(String message) {
	            super(String.format("The provided VCF file is malformed: %s", message));
	        }

	        public MalformedVCF(String message, int lineNo) {
	            super(String.format("The provided VCF file is malformed at approximately line number %d: %s", lineNo, message));
	        }
    }

    public static class MalformedBAM extends UserException {
        /**
		 * 
		 */
		private static final long serialVersionUID = 1476767881468371580L;

		public MalformedBAM(SAMRecord read, String message) {
            this(read.getFileSource() != null ? read.getFileSource().getReader().toString() : "(none)", message);
        }

        public MalformedBAM(File file, String message) {
            this(file.toString(), message);
        }

        public MalformedBAM(String source, String message) {
            super(String.format("SAM/BAM file %s is malformed: %s", source, message));
        }
    }
    
    public static class ReadMissingReadGroup extends MalformedBAM {
        /**
		 * 
		 */
		private static final long serialVersionUID = 7071617091686069309L;

		public ReadMissingReadGroup(SAMRecord read) {
            super(read, String.format("Read %s is either missing the read group or its read group is not defined in the BAM header, both of which are required by the GATK.  Please use http://www.broadinstitute.org/gsa/wiki/index.php/ReplaceReadGroups to fix this problem", read.getReadName()));
        }
    }
    
    public static class CommandLineException extends UserException {
        /**
		 * 
		 */
		private static final long serialVersionUID = 5637179533604388076L;

		public CommandLineException(String message) {
            super(String.format("Invalid command line: %s", message));
        }
    }
 
    public static class BadArgumentValueException extends CommandLineException {
        /**
		 * 
		 */
		private static final long serialVersionUID = -2902450578967979989L;

		public BadArgumentValueException(String arg, String message) {
            super(String.format("Argument %s has a bad value: %s", arg, message));
        }
		
		public BadArgumentValueException(int value,int min,int max){
			super(String.format("Argument must in (%d,%d),but user is set to %d", min,max,value));
		}
		
		public BadArgumentValueException(String arg, String value, String message){
            super(String.format("Argument %s has a bad value: %s. %s", arg, value,message));
        }
    }
    
    public static class CouldNotCreateOutputFile extends UserException {
        /**
		 * 
		 */
		private static final long serialVersionUID = -6518201215435711936L;

		public CouldNotCreateOutputFile(File file, String message, Exception e) {
            super(String.format("Couldn't write file %s because %s with exception %s", file.getAbsolutePath(), message, getMessage(e)));
        }

        public CouldNotCreateOutputFile(File file, String message) {
            super(String.format("Couldn't write file %s because %s", file.getAbsolutePath(), message));
        }

        public CouldNotCreateOutputFile(String file, String message) {
            super(String.format("Couldn't write file %s because %s", file, message));
        }

        
        public CouldNotCreateOutputFile(String filename, String message, Exception e) {
            super(String.format("Couldn't write file %s because %s with exception %s", filename, message, getMessage(e)));
        }

        public CouldNotCreateOutputFile(File file, Exception e) {
            super(String.format("Couldn't write file %s because exception %s", file.getAbsolutePath(), getMessage(e)));
        }

        public CouldNotCreateOutputFile(String message, Exception e) {
            super(message, e);
        }
    }

    
    public static class CouldNotReadInputFile extends UserException {
        /**
		 * 
		 */
		private static final long serialVersionUID = -7414815689723053752L;

		public CouldNotReadInputFile(String message, Exception e) {
            super(String.format("Couldn't read file because %s caused by %s", message, getMessage(e)));
        }

        public CouldNotReadInputFile(File file) {
            super(String.format("Couldn't read file %s", file.getAbsolutePath()));
        }

        public CouldNotReadInputFile(File file, String message) {
            super(String.format("Couldn't read file %s because %s", file.getAbsolutePath(), message));
        }
        public CouldNotReadInputFile(String file, String message) {
            super(String.format("Couldn't read file %s because %s", file, message));
        }
        public CouldNotReadInputFile(File file, String message, Exception e) {
            super(String.format("Couldn't read file %s because %s with exception %s", file.getAbsolutePath(), message, getMessage(e)));
        }

        public CouldNotReadInputFile(File file, Exception e) {
            this(file, getMessage(e));
        }

        public CouldNotReadInputFile(String message) {
            super(message);
        }
    }


    public static class ReadAttributeTypeMismatch extends UserException {
        private static final long serialVersionUID = 0L;

        public ReadAttributeTypeMismatch( final String attributeName, final String targetType ) {
            super(String.format("Attribute %s not of (or convertible to) type %s", attributeName, targetType));
        }

        public ReadAttributeTypeMismatch( final String attributeName, final String targetType, final Throwable throwable ) {
            super(String.format("Attribute %s not of (or convertible to) type %s", attributeName, targetType), throwable);
        }

        public ReadAttributeTypeMismatch( final SAMRecord read, final String attributeName, final String targetType) {
            super(String.format("In read %s @ %s attribute %s not of (or convertible to) type %s", read.getReadName(), "" + read.getContig() + ":" + read.getStart(), attributeName, targetType));
        }

        public ReadAttributeTypeMismatch( final SAMRecord read, final String attributeName, final String targetType, final Object value) {
            super(String.format("In read %s @ %s attribute %s not of (or convertible to) type %s: class is '%s' and value is '%s'", read.getReadName(), "" + read.getContig() + ":" + read.getStart(), attributeName, targetType, value == null ? "null" : value.getClass(), "" + value));
        }

        public ReadAttributeTypeMismatch( final SAMRecord read, final String attributeName, final String targetType, final Throwable ex) {
            super(String.format("In read %s @ %s attribute %s not of (or convertible to) type %s", read.getReadName(), "" + read.getContig() + ":" + read.getStart(), attributeName, targetType, ex));
        }
    }
}
