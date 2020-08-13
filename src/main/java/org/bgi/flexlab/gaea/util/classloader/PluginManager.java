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

package org.bgi.flexlab.gaea.util.classloader;



import org.bgi.flexlab.gaea.data.exception.DynamicClassResolutionException;
import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.tools.genotyer.annotator.*;
import org.bgi.flexlab.gaea.tools.genotyer.annotator.interfaces.*;

import java.lang.reflect.Constructor;
import java.lang.reflect.Method;
import java.net.URL;
import java.net.URLClassLoader;
import java.util.*;

public class PluginManager<PluginType> {

    /**
     * Defines the category of plugin defined by the subclass.
     */
    protected final String pluginCategory;

    /**
     * Define common strings to trim off the end of the name.
     */
    protected final String pluginSuffix;
    
    /**
     * Plugins stored based on their name.
     */
    private final SortedMap<String, Class<? extends PluginType>> pluginsByName;

    private final List<Class<? extends PluginType>> plugins;
    private final List<Class<? extends PluginType>> interfaces;
    private static List<Class> allClass;
    static{
    	allClass=new ArrayList<Class>();
        allClass.add(AlleleBalance.class);
        allClass.add(AlleleBalanceBySample.class);
    	allClass.add(BaseCounts.class);
    	allClass.add(BaseQualityRankSumTest.class);
    	allClass.add(ChromosomeCounts.class);
    	allClass.add(ClippingRankSumTest.class);
    	allClass.add(DepthOfCoverage.class);
    	allClass.add(DepthPerAlleleBySample.class);
    	allClass.add(FisherStrand.class);
    	allClass.add(GCContent.class);
    	allClass.add(HaplotypeScore.class);
    	allClass.add(HardyWeinberg.class);
    	allClass.add(HomopolymerRun.class);
    	allClass.add(InbreedingCoeff.class);
    	allClass.add(IndelType.class);
    	allClass.add(LowMQ.class);
    	allClass.add(MappingQualityRankSumTest.class);
    	allClass.add(MappingQualityZero.class);
    	allClass.add(MappingQualityZeroBySample.class);
    	allClass.add(MappingQualityZeroFraction.class);
    	//allClass.add(MVLikelihoodRatio.class);
    	//allClass.add(NBaseCount.class);
    	allClass.add(QualByDepth.class);
    	allClass.add(RankSumTest.class);
    	allClass.add(ReadPosRankSumTest.class);
    	allClass.add(RMSMappingQuality.class);
    	allClass.add(SampleList.class);
    	//allClass.add(SnpEff.class);
    	allClass.add(SpanningDeletions.class);
    	allClass.add(TandemRepeatAnnotator.class);
    	//allClass.add(TechnologyComposition.class);
    	//allClass.add(TransmissionDisequilibriumTest.class);
    	allClass.add(DBsnp.class);
    	allClass.add(VariantAnnotatorEngine.class);
    	allClass.add(ActiveRegionBasedAnnotation.class);
    	allClass.add(AnnotationInterfaceManager.class);
    	allClass.add(AnnotationType.class);
    	allClass.add(ExperimentalAnnotation.class);
    	allClass.add(GenotypeAnnotation.class);
    	allClass.add(InfoFieldAnnotation.class);
    	allClass.add(RodRequiringAnnotation.class);
    	allClass.add(StandardAnnotation.class);//
    	allClass.add(VariantAnnotatorAnnotation.class);
    	allClass.add(WorkInProgressAnnotation.class);
    }
    /**
     * Create a new plugin manager.
     * @param pluginType Core type for a plugin.
     */
    public PluginManager(Class pluginType) {
        this(pluginType, pluginType.getSimpleName().toLowerCase(), pluginType.getSimpleName(), null);
    }

    /**
     * Create a new plugin manager.
     * @param pluginType Core type for a plugin.
     * @param classpath Custom class path to search for classes.
     */
    public PluginManager(Class pluginType, List<URL> classpath) {
        this(pluginType, pluginType.getSimpleName().toLowerCase(), pluginType.getSimpleName(), classpath);
    }

    /**
     * Create a new plugin manager.
     * @param pluginType Core type for a plugin.
     * @param pluginCategory Provides a category name to the plugin.  Must not be null.
     * @param pluginSuffix Provides a suffix that will be trimmed off when converting to a plugin name.  Can be null.
     */
    public PluginManager(Class pluginType, String pluginCategory, String pluginSuffix) {
        this(pluginType, pluginCategory, pluginSuffix, null);
    }

    /**
     * Create a new plugin manager.
     * @param pluginType Core type for a plugin.
     * @param pluginCategory Provides a category name to the plugin.  Must not be null.
     * @param pluginSuffix Provides a suffix that will be trimmed off when converting to a plugin name.  Can be null.
     * @param classpath Custom class path to search for classes.
     */
    public PluginManager(Class pluginType, String pluginCategory, String pluginSuffix, List<URL> classpath) {
        this.pluginCategory = pluginCategory;
        this.pluginSuffix = pluginSuffix;
        this.plugins = new ArrayList<Class<? extends PluginType>>();
        this.interfaces = new ArrayList<Class<? extends PluginType>>();

        

        // Load all classes types filtering them by concrete.
        //@SuppressWarnings("unchecked")
       // Set<Class<? extends PluginType>> allTypes = reflections.getSubTypesOf(pluginType);
        Set<Class<? extends PluginType>> allTypes=loadClass(pluginType);
        //System.out.println("allType:"+allTypes.size());
        //for(Class c:allTypes)
        //	System.out.println(c.toString());
        for( Class<? extends PluginType> type: allTypes ) {
            // The plugin manager does not support anonymous classes; to be a plugin, a class must have a name.
            if(JVMUtils.isAnonymous(type))
                continue;

            if( JVMUtils.isConcrete(type) )
                plugins.add(type);
            else
                interfaces.add(type);
        }

        pluginsByName = new TreeMap<String, Class<? extends PluginType>>();
        for (Class<? extends PluginType> pluginClass : plugins) {
            String pluginName = getName(pluginClass);
            pluginsByName.put(pluginName, pluginClass);
        }

        // sort the plugins so the order of elements is deterministic
        sortPlugins(plugins);
        sortPlugins(interfaces);
    }

    /**
     * Sorts, in place, the list of plugins according to getName() on each element
     *
     * @param unsortedPlugins unsorted plugins
     */
    private void sortPlugins(final List<Class<? extends PluginType>> unsortedPlugins) {
        Collections.sort(unsortedPlugins, new ComparePluginsByName());
    }

    private final class ComparePluginsByName implements Comparator<Class<? extends PluginType>> {
        @Override
        public int compare(final Class<? extends PluginType> aClass, final Class<? extends PluginType> aClass1) {
            String pluginName1 = getName(aClass);
            String pluginName2 = getName(aClass1);
            return pluginName1.compareTo(pluginName2);
        }
    }

    /**
     * Adds the URL to the system class loader classpath using reflection.
     * HACK: Uses reflection to modify the class path, and assumes loader is a URLClassLoader.
     * @param urls URLs to add to the system class loader classpath.
     */
    private static void addClasspath(List<URL> urls) {
      Collection<URL> existing = JVMUtils.getClasspathURLs();
      for (URL url : urls) {
          if (existing.contains(url))
            continue;
          try {
              Method method = URLClassLoader.class.getDeclaredMethod("addURL", URL.class);
              if (!method.isAccessible())
                  method.setAccessible(true);
              method.invoke(ClassLoader.getSystemClassLoader(), url);
          } catch (Exception e) {
              throw new UserException("Error adding url to the current classloader.", e);
          }
      }
    }
    
    public Map<String, Class<? extends PluginType>> getPluginsByName() {
        return Collections.unmodifiableMap(pluginsByName);
    }

    /**
     * Does a plugin with the given name exist?
     *
     * @param pluginName Name of the plugin for which to search.
     * @return True if the plugin exists, false otherwise.
     */
    public boolean exists(String pluginName) {
        return pluginsByName.containsKey(pluginName);
    }

    /**
     * Does a plugin with the given name exist?
     *
     * @param plugin Name of the plugin for which to search.
     * @return True if the plugin exists, false otherwise.
     */
    public boolean exists(Class<? extends PluginType> plugin) {
        return pluginsByName.containsValue(plugin);
    }

    /**
     * Returns the plugin classes
     * @return the plugin classes
     */
    public List<Class<? extends PluginType>> getPlugins() {
        return plugins;
    }

    /**
     * Returns the interface classes
     * @return the interface classes
     */
    public List<Class<? extends PluginType>> getInterfaces() {
        return interfaces;
    }

    /**
     * Returns the plugin classes implementing interface or base clase
     * @param type type of interface or base class
     * @return the plugin classes implementing interface or base class
     */
    public List<Class<? extends PluginType>> getPluginsImplementing(Class<?> type) {
        List<Class<? extends PluginType>> implementing = new ArrayList<Class<? extends PluginType>>();
        for (Class<? extends PluginType> plugin: getPlugins())
            if (type.isAssignableFrom(plugin))
                implementing.add(plugin);
        return implementing;
    }



    /**
     * Gets a plugin with the given name
     *
     * @param pluginName Name of the plugin to retrieve.
     * @return The plugin object if found; null otherwise.
     */
    public PluginType createByName(String pluginName) {
        Class<? extends PluginType> plugin = pluginsByName.get(pluginName);
        if( plugin == null ) {
            String errorMessage = formatErrorMessage(pluginCategory,pluginName);
        }
        try {
            return plugin.newInstance();
        } catch (Exception e) {
            throw new DynamicClassResolutionException(plugin, e);
        }
    }

    /**
     * create a plugin with the given type
     *
     * @param pluginType type of the plugin to create.
     * @return The plugin object if created; null otherwise.
     */
    public PluginType createByType(Class<? extends PluginType> pluginType) {
        try {
            Constructor<? extends PluginType> noArgsConstructor = pluginType.getDeclaredConstructor((Class[])null);
            noArgsConstructor.setAccessible(true);
            return noArgsConstructor.newInstance();
        } catch (Exception e) {
            throw new DynamicClassResolutionException(pluginType, e);
        }
    }

    /**
     * Returns concrete instances of the plugins
     * @return concrete instances of the plugins
     */
    public List<PluginType> createAllTypes() {
        List<PluginType> instances = new ArrayList<PluginType>();
        for ( Class<? extends PluginType> c : getPlugins() ) {
            instances.add(createByType(c));
        }
        return instances;
    }

    /**
     * Create a name for this type of plugin.
     *
     * @param pluginType The type of plugin.
     * @return A name for this type of plugin.
     */
    public String getName(Class pluginType) {
        String pluginName = "";

        if (pluginName.length() == 0) {
            pluginName = pluginType.getSimpleName();
            if (pluginSuffix != null && pluginName.endsWith(pluginSuffix))
                pluginName = pluginName.substring(0, pluginName.lastIndexOf(pluginSuffix));
        }

        return pluginName;
    }
    /**
     * load class 
     */
    protected Set<Class<? extends PluginType>> loadClass(Class<PluginType> pluginType)
    {
    	
    	Set<Class<? extends PluginType>> typeClass=new HashSet<Class<? extends PluginType>>();
    	for(Class c:allClass)
    		if(pluginType.isAssignableFrom(c))
    		{
    			if(c.getName().equals(pluginType.getName()))	
    				continue;
    			typeClass.add(c);
    		}
    	return typeClass;
    }
    /**
     * Generate the error message for the plugin manager. The message is allowed to depend on the class.
     * @param pluginCategory - string, the category of the plugin (e.g. read filter)
     * @param pluginName - string, what we were trying to match (but failed to)
     * @return error message text describing the error
     */
    protected String formatErrorMessage(String pluginCategory, String pluginName ) {
        return String.format("Could not find %s with name: %s", pluginCategory,pluginName);
    }
    
    public static void main(String args[]){
        //Class<InfoFieldAnnotation> info=InfoFieldAnnotation.class;
        //Class<SampleList> sample=SampleList.class;
       // if(info.isAssignableFrom( sample))
        //	System.out.println("==================");
    	//if(InfoFieldAnnotation)
    	PluginManager<InfoFieldAnnotation> infoFieldAnnotationPluginManager = new PluginManager<InfoFieldAnnotation>(InfoFieldAnnotation.class);
   	 //PluginManager<GenotypeAnnotation> genotypeAnnotationPluginManager = new PluginManager<GenotypeAnnotation>(GenotypeAnnotation.class);
    	//PluginManager<AnnotationType> annotationTypePluginManager = new PluginManager<AnnotationType>(AnnotationType.class);
     }
    
}
