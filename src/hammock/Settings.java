/* Singleton.
 * Class loads settings saved in configuration files to be able to serve them
 * to methods later. 
 */
package hammock;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Properties;

/**
 *
 * @author Adam Krejci
 */
public class Settings {

    private static final Settings instance = new Settings();
    private final char separatorChar = File.separatorChar; //system separator
    private final String jarParentDir;
    private final String propertiesFile;

    private String clustalCommand = null;
    private String hmmbuildCommand = null;
    private String hmmsearchCommand = null;
    private String hmmalignCommand = null;
    private String hhmakeCommand = null;
    private String hhsearchCommand = null;

    private List<String> clustalParameters = null;
    private List<String> hmmbuildParameters = null;
    private List<String> hmmsearchParameters = null;
    private List<String> hmmalignParameters = null;
    private List<String> hhmakeParameters = null;
    private List<String> hhsearchParameters = null;


    private String tempPath = null;
    private static String tempName;

    private String fastaDirectory = null;
    private String msaDirectory = null;
    private String hmmDirectory = null;
    private String hhDirectory = null;
    private String hmmsearchOutDirectory = null;
    private String hhsearchOutDirectory = null;
    private String fastaDatabaseFile = null;

    public static Settings getInstance() {
        return instance;
    }

    /**
     * Constructor loads settings file. If settings file is not found, a warning
     * is printed in stderr and all settings attributes remain null.
     */
    private Settings() {
        jarParentDir = (new File(Hammock.class.getProtectionDomain().getCodeSource()
                .getLocation().getPath()).getParentFile().getParentFile().getPath());
        propertiesFile = jarParentDir + separatorChar + "settings" + separatorChar + "settings.prop";
        tempName = "Hammock_temp_" + System.currentTimeMillis();

        try {
            Properties properties = loadProperties(propertiesFile);
            clustalCommand = properties.getProperty("clustalOmegaCommand",
                    jarParentDir + separatorChar + "clustal-omega-1.2.0" + separatorChar + "clustalO-64bit");
            hmmbuildCommand = properties.getProperty("hmmbuildCommand",
                    jarParentDir + separatorChar + "hmmer-3.1b1" + separatorChar + "src" + separatorChar + "hmmbuild");
            hmmsearchCommand = properties.getProperty("hmmsearchCommand",
                    jarParentDir + separatorChar + "hmmer-3.1b1" + separatorChar + "src" + separatorChar + "hmmsearch");
            hmmalignCommand = properties.getProperty("hmmalignCommand", 
                    jarParentDir + separatorChar + "hmmer-3.1b1" + separatorChar + "src" + separatorChar + "hmmalign");
            hhmakeCommand = properties.getProperty("hhmakeCommand",
                    jarParentDir + separatorChar + "hhsuite-2.0.16" + separatorChar + "bin" + separatorChar + "hhmake");
            hhsearchCommand = properties.getProperty("hhsearchCommand",
                    jarParentDir + separatorChar + "hhsuite-2.0.16" + separatorChar + "bin" + separatorChar + "hhsearch");

            clustalParameters = loadListOfParameters(properties, "clustalOmegaParameters");
            hmmbuildParameters = loadListOfParameters(properties, "hmmbuildParameters");
            hmmsearchParameters = loadListOfParameters(properties, "hmmsearchParameters");
            hmmalignParameters = loadListOfParameters(properties, "hmmalignParameters");
            hhmakeParameters = loadListOfParameters(properties, "hhmakeParameters");
            hhsearchParameters = loadListOfParameters(properties, "hhsearchParameters");

            List<String> directoriesToCreate = new ArrayList<>();

            String tempDirectory = properties.getProperty("tempDirectory", "/tmp");
            tempPath = tempDirectory + separatorChar + tempName;
            directoriesToCreate.add(tempPath);
            fastaDirectory = properties.getProperty("fastaDirectory", tempPath + separatorChar + "fasta" + separatorChar);
            directoriesToCreate.add(fastaDirectory);
            msaDirectory = properties.getProperty("msaDirectory", tempPath + separatorChar + "msa" + separatorChar);
            directoriesToCreate.add(msaDirectory);
            hmmDirectory = properties.getProperty("hmmDirectory", tempPath + separatorChar + "hmm" + separatorChar);
            directoriesToCreate.add(hmmDirectory);
            hhDirectory = properties.getProperty("hhDirectory", tempPath + separatorChar + "hh" + separatorChar);
            directoriesToCreate.add(hhDirectory);
            hmmsearchOutDirectory = properties.getProperty("hmmsearchOutDirectory", tempPath + separatorChar + "hmmsearchOut" + separatorChar);
            directoriesToCreate.add(hmmsearchOutDirectory);
            hhsearchOutDirectory = properties.getProperty("hhsearchOutDirectory", tempPath + separatorChar + "hhsearchOut" + separatorChar);
            directoriesToCreate.add(hhsearchOutDirectory);
            fastaDatabaseFile = properties.getProperty("fastaDatabaseFile", tempPath + separatorChar + "fastaDatabase.fa");

            for (String path : directoriesToCreate) {
                File file = new File(path);
                if (!file.exists()) {
                    file.mkdir();
                }
            }

        } catch (IOException e) {
            System.err.println("Warning: Settings file not loaded."
                    + " Any attempt to use clustal, hmmer or hhsuite will fail. "
                    + "Properties file should be in: " + propertiesFile + ". "
                    + "Error message: " + e.getMessage());

        }
    }

    /**
     * Loads property for specified key and parses its contents as
     * whitespace-delimited list of program parameters. Returns resulting list.
     *
     * @param properties Properties object from which property will be loaded
     * @param propertyKey Key under which property to be parsed will be selected
     * @return List of parameters from property specified
     */
    private List<String> loadListOfParameters(Properties properties, String propertyKey) {
        String paramString = properties.getProperty(propertyKey);
        List<String> result = null;
        if (paramString != null) {
            result = new ArrayList<>();
            for (String param : paramString.split("\\s+")) {//split at any number of any whitespace characters
                if (param.length() > 0) {
                    result.add(param);
                }
            }
        }
        return result;
    }

    /**
     * Loads properties file into a Properties object.
     *
     * @param fileName Path to file with properties.
     * @return Properties object containing properties from file specified by
     * fileName.
     * @throws IOException
     */
    private Properties loadProperties(String fileName) throws IOException {
        try (FileReader fReader = new FileReader(new File(fileName))) {
            Properties properties = new Properties();
            properties.load(fReader);
            return properties;
        } catch (IOException e) {
            throw new IOException(e);
        }
    }

    public char getSeparatorChar() {
        return separatorChar;
    }

    public String getClustalCommand() {
        return clustalCommand;
    }

    public String getHmmbuildCommand() {
        return hmmbuildCommand;
    }

    public String getHmmsearchCommand() {
        return hmmsearchCommand;
    }

    public String getHmmalignCommand() {
        return hmmalignCommand;
    }
    

    public String getHhmakeCommand() {
        return hhmakeCommand;
    }


    public String getHhsearchCommand() {
        return hhsearchCommand;
    }

    public List<String> getClustalParameters() {
        return clustalParameters;
    }

    public List<String> getHmmbuildParameters() {
        return hmmbuildParameters;
    }

    public List<String> getHmmsearchParameters() {
        return hmmsearchParameters;
    }

    public List<String> getHmmalignParameters() {
        return hmmalignParameters;
    }
    
    public List<String> getHhmakeParameters() {
        return hhmakeParameters;
    }


    public List<String> getHhsearchParameters() {
        return hhsearchParameters;
    }


    public String getFastaDirectory() {
        return fastaDirectory;
    }

    public String getMsaDirectory() {
        return msaDirectory;
    }

    public String getHmmDirectory() {
        return hmmDirectory;
    }

    public String getHmmsearchOutDirectory() {
        return hmmsearchOutDirectory;
    }


    public String getHhsearchOutDirectory() {
        return hhsearchOutDirectory;
    }

    public String getFastaDatabaseFile() {
        return fastaDatabaseFile;
    }

    public String getHhDirectory() {
        return hhDirectory;
    }

    public String getTempDirectory() {
        return tempPath;
    }

}
