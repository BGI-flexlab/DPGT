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
package org.bgi.flexlab.gaea.tools.genotyer.annotator.interfaces;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.TreeSet;

import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.util.classloader.PluginManager;


public class AnnotationInterfaceManager {
    private static PluginManager<InfoFieldAnnotation> infoFieldAnnotationPluginManager = new PluginManager<InfoFieldAnnotation>(InfoFieldAnnotation.class);
    private static PluginManager<GenotypeAnnotation> genotypeAnnotationPluginManager = new PluginManager<GenotypeAnnotation>(GenotypeAnnotation.class);
    private static PluginManager<AnnotationType> annotationTypePluginManager = new PluginManager<AnnotationType>(AnnotationType.class);

    public static List<InfoFieldAnnotation> createAllInfoFieldAnnotations() {
    	return infoFieldAnnotationPluginManager.createAllTypes();
    }

    public static List<GenotypeAnnotation> createAllGenotypeAnnotations() {
    	//System.out.println(genotypeAnnotationPluginManager.createAllTypes().size());
    	//for(GenotypeAnnotation  info:genotypeAnnotationPluginManager.createAllTypes())
    	//	System.out.println(info.getClass());
    	return genotypeAnnotationPluginManager.createAllTypes();
    }

    public static void validateAnnotations(List<String> annotationGroupsToUse, List<String> annotationsToUse) {
        HashMap<String, Class> classMap = new HashMap<String, Class>();
        for ( Class c : infoFieldAnnotationPluginManager.getPlugins() )
        {
        	//System.out.println(c.getSimpleName()+"\t"+c.toString());

        	classMap.put(c.getSimpleName(), c);
        }
        for ( Class c : genotypeAnnotationPluginManager.getPlugins() )
        {
        	//System.out.println(c.getSimpleName()+"\t"+c.toString());

        	classMap.put(c.getSimpleName(), c);
        }
        for ( Class c : annotationTypePluginManager.getInterfaces() )
        {
        	//System.out.println(c.getSimpleName()+"\t"+c.toString());

        	classMap.put(c.getSimpleName(), c);
        }

        if ( annotationGroupsToUse.size() != 1 || !"none".equals(annotationGroupsToUse.get(0)) ) {
            for ( String group : annotationGroupsToUse ) {
                Class interfaceClass = classMap.get(group);
                if ( interfaceClass == null )
                    interfaceClass = classMap.get(group + "Annotation");
                if ( interfaceClass == null )
                    throw new UserException.BadArgumentValueException("group", "Class " + group + " is not found; please check that you have specified the class name correctly");
            }
        }

        // validate the specific classes provided
        for ( String annotation : annotationsToUse ) {
            Class annotationClass = classMap.get(annotation);
            if ( annotationClass == null )
                annotationClass = classMap.get(annotation + "Annotation");
            if ( annotationClass == null )
                throw new UserException.BadArgumentValueException("annotation", "Class " + annotation + " is not found; please check that you have specified the class name correctly");
        }
    }

    public static List<InfoFieldAnnotation> createInfoFieldAnnotations(List<String> annotationGroupsToUse, List<String> annotationsToUse) {
       // System.out.println(createAnnotations(infoFieldAnnotationPluginManager, annotationGroupsToUse, annotationsToUse).size());
    	return createAnnotations(infoFieldAnnotationPluginManager, annotationGroupsToUse, annotationsToUse);
    }

    public static List<GenotypeAnnotation> createGenotypeAnnotations(List<String> annotationGroupsToUse, List<String> annotationsToUse) {
        return createAnnotations(genotypeAnnotationPluginManager, annotationGroupsToUse, annotationsToUse);
    }

    private static <T> List<T> createAnnotations(PluginManager<T> pluginManager, List<String> annotationGroupsToUse, List<String> annotationsToUse) {
        // get the instances
        List<T> annotations = new ArrayList<T>();

        // get the classes from the provided groups (interfaces)
        // create a map for all annotation classes which implement our top-level interfaces
        HashMap<String, Class> classMap = new HashMap<String, Class>();
        for ( Class c : pluginManager.getPlugins() )
        	classMap.put(c.getSimpleName(), c);
        for ( Class c : annotationTypePluginManager.getInterfaces() )
        	classMap.put(c.getSimpleName(), c);
 
        // use a TreeSet so that classes are returned deterministically (the plugin manager apparently isn't deterministic)
        TreeSet<Class> classes = new TreeSet<Class>(new Comparator<Class>() {
            public int compare(Class o1, Class o2) {
                return o1.getSimpleName().compareTo(o2.getSimpleName());
            }
        });

        if ( annotationGroupsToUse.size() != 1 || !"none".equals(annotationGroupsToUse.get(0)) ) {
            for ( String group : annotationGroupsToUse ) {
                Class interfaceClass = classMap.get(group);
                if ( interfaceClass == null )
                    interfaceClass = classMap.get(group + "Annotation");
                if ( interfaceClass != null )
                {
                	//for(Class c:pluginManager.getPluginsImplementing(interfaceClass))
                	//	System.out.println(c.toString());
                	classes.addAll(pluginManager.getPluginsImplementing(interfaceClass));
                }
            }
        }

        // get the specific classes provided
        for ( String annotation : annotationsToUse ) {
            Class annotationClass = classMap.get(annotation);
            if ( annotationClass == null )
                annotationClass = classMap.get(annotation + "Annotation");
            if ( annotationClass != null )
                classes.add(annotationClass);
        }

        // note that technically an annotation can work on both the INFO and FORMAT fields
        for ( Class c : classes )
            annotations.add((T)pluginManager.createByType(c));
        //System.out.println(annotations.size());
        return annotations;
    }
    
    public static void main(String args[]){
    	List<String> toUse=new ArrayList<String>();
    	toUse.add("Standard");
    	//AnnotationInterfaceManager manager=new AnnotationInterfaceManager();
    	AnnotationInterfaceManager.createInfoFieldAnnotations(toUse,new ArrayList<String>());
    }
    
}
