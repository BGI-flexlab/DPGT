package org.bgi.flexlab.gaea.util;

import java.io.File;
import java.lang.reflect.Constructor;
import java.lang.reflect.Method;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLClassLoader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.bgi.flexlab.gaea.data.exception.DynamicClassResolutionException;
import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.util.classloader.JVMUtils;
import org.reflections.Reflections;
import org.reflections.scanners.SubTypesScanner;
import org.reflections.util.ConfigurationBuilder;

public class PluginManager<PluginType> {
    /**
     * A reference into our introspection utility.
     */
    private static final Reflections defaultReflections;

    static {
        // turn off logging in the reflections library - they talk too much
        Reflections.log = null;

        Set<URL> classPathUrls = new LinkedHashSet<URL>();

        URL cwd;
        try {
            cwd = new File("").getAbsoluteFile().toURI().toURL();
        } catch (MalformedURLException e) {
            throw new RuntimeException(e);
        }

        // NOTE: Reflections also scans directories for classes.
        // Meanwhile some of the jar MANIFEST.MF Bundle-ClassPath properties contain "."
        // Do NOT let reflections scan the CWD where it often picks up test classes when
        // they weren't explicitly in the classpath, for example the UninstantiableWalker
        for (URL url: JVMUtils.getClasspathURLs())
            if (!url.equals(cwd))
                classPathUrls.add(url);

        defaultReflections = new Reflections( new ConfigurationBuilder()
            .setUrls(classPathUrls)
            .setScanners(new SubTypesScanner()));
    }

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

        Reflections reflections;
        if (classpath == null) {
            reflections = defaultReflections;
        } else {
            addClasspath(classpath);
            reflections = new Reflections( new ConfigurationBuilder()
                .setUrls(classpath)
                .setScanners(new SubTypesScanner()));
        }

        // Load all classes types filtering them by concrete.
        @SuppressWarnings("unchecked")
        Set<Class<? extends PluginType>> allTypes = reflections.getSubTypesOf(pluginType);
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
              throw new RuntimeException("Error adding url to the current classloader.");
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
            throw createMalformedArgumentException(errorMessage);
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
        Logger logger = Logger.getLogger(PluginManager.class);
        logger.setLevel(Level.ERROR);
        try {
            Constructor<? extends PluginType> noArgsConstructor = pluginType.getDeclaredConstructor((Class[])null);
            noArgsConstructor.setAccessible(true);
            return noArgsConstructor.newInstance();
        } catch (Exception e) {
            logger.error("Couldn't initialize the plugin. Typically this is because of wrong global class variable initializations.");
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
     * Generate the error message for the plugin manager. The message is allowed to depend on the class.
     * @param pluginCategory - string, the category of the plugin (e.g. read filter)
     * @param pluginName - string, what we were trying to match (but failed to)
     * @return error message text describing the error
     */
    protected String formatErrorMessage(String pluginCategory, String pluginName ) {
        return String.format("Could not find %s with name: %s", pluginCategory,pluginName);
    }

    /**
     * Creates a UserException with the appropriate message for this instance.
     * @param errorMessage formatted error message from formatErrorMessage().
     * @return A UserException with the error message.
     */
    protected UserException createMalformedArgumentException(final String errorMessage) {
        throw new UserException(errorMessage);
    }
}

