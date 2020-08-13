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
 *
 * This file incorporates work covered by the following copyright and 
 * Permission notices:
 *
 * Copyright (c) 2009-2012 The Broad Institute
 *  
 *     Permission is hereby granted, free of charge, to any person
 *     obtaining a copy of this software and associated documentation
 *     files (the "Software"), to deal in the Software without
 *     restriction, including without limitation the rights to use,
 *     copy, modify, merge, publish, distribute, sublicense, and/or sell
 *     copies of the Software, and to permit persons to whom the
 *     Software is furnished to do so, subject to the following
 *     conditions:
 *  
 *     The above copyright notice and this permission notice shall be
 *     included in all copies or substantial portions of the Software.
 *  
 *     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 *     OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *     NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 *     WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 *     FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 *     OTHER DEALINGS IN THE SOFTWARE.
 *******************************************************************************/
package org.bgi.flexlab.gaea.tools.vcfqualitycontrol.variantrecalibratioin.traindata;

import htsjdk.samtools.util.RuntimeEOFException;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.bgi.flexlab.gaea.data.structure.header.VCFConstants;

import java.sql.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class KnownDatabase {
	private static final String DBDRIVER = "com.mysql.jdbc.Driver";
    private static final String DBURL = "jdbc:mysql://10.1.0.60:3306/KnownSNP";
	//private static final String DBURL = "jdbc:mysql://172.30.4.57:3306/KnownSNP";
    private static final String DBUSER = "gaea";
    private static final String DBPASSWORD = "gaea";
    private static Connection conn;
    private static Statement statement;
    private static int lineNo;

    public KnownDatabase() throws ClassNotFoundException {
    	Class.forName(DBDRIVER);
    	try {
			conn = DriverManager.getConnection(DBURL, DBUSER, DBPASSWORD);
			if(!conn.isClosed()) 
	             System.out.println("Succeeded connecting to the Database!");
	    	 statement = conn.createStatement();
		} catch (SQLException e) {
			throw new RuntimeException(e);
		}
    	 
    	 lineNo = 0;
    }
    public ResultSet getSourceId(String refType, String name) throws SQLException {
    	StringBuilder sql = new StringBuilder();
    	sql.append("select * from source where ");
    	if(refType != null) {
    		sql.append("ref_type='");
    		sql.append(refType);
    		sql.append("'");
    	}
    	if(name != null) {
    		if(refType != null)
    			sql.append(" && ");
    		sql.append("name='");
    		sql.append(name);
    		sql.append("'");
    	}
    	sql.append(";");
    	return statement.executeQuery(sql.toString());
    }

    
    public ArrayList<VariantContext> getTrainData(int resourceId, String chrName, int start, int end) {
    	ArrayList<VariantContext> vcl = new ArrayList<VariantContext>();
        ResultSet rs;
		try {
			StringBuilder sql = new StringBuilder();
			sql.append("select * from variation where ");
			sql.append("source_id=");
			sql.append(resourceId);
			sql.append(" && chr='");
			sql.append(chrName);
			sql.append("' && pos=");
			sql.append(start);
			sql.append(";");
			rs = statement.executeQuery(sql.toString());
			while(rs.next()) {
				VariantContextBuilder vcb = new VariantContextBuilder();
		       	buildVariationContext(vcb, rs);
		       	VariantContext vc = vcb.make();
		       	if(vc.getEnd() <= end)
		       		vcl.add(vc);
			 }
		} catch (SQLException e) {
			throw new RuntimeEOFException(e.toString());
		}     
    	
    	return vcl;
    }
    
    private void buildVariationContext(VariantContextBuilder vcb, ResultSet rs) throws SQLException {
    	lineNo++;
    	vcb.chr(rs.getString(2));
    	vcb.start(rs.getInt(3));
    	vcb.id(rs.getString(4));
    	
    	Allele ref = Allele.create(rs.getString(5), true);
    	Allele alt = Allele.create(rs.getString(6), false);
    	final List<Allele> alleles = new ArrayList<>();
    	Collections.addAll(alleles, ref, alt);
    	vcb.alleles(alleles);
    	
    	vcb.stop(rs.getInt(3) + rs.getString(5).length() - 1);
    	String filterString = rs.getString(7);
    	if(!filterString.equals(VCFConstants.UNFILTERED) && !filterString.equals(VCFConstants.PASSES_FILTERS_v4))
    		vcb.filter(filterString);
    }
}
