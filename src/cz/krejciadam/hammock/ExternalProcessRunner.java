/* 
 * Class executes programs as external processes.
 */

package cz.krejciadam.hammock;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 *
 * @author Adam Krejci
 */
public class ExternalProcessRunner {
    
    /**
     * Runs an external process and waits for its completion. Writes stdout
     * and stderr of that process to printStreams specified. 
     * @param programCommand A command to be called as an external process
     * @param arguments Program arguments in String form. Each entry (separated
     * by space when calling from command line) should be a single String in this
     * list. E.g. "-nthreads 8" should be two separate Strings without any gaps
     * in parameter arguments.
     * @param out a PrintStream object to print stdout to
     * @param errOut a PrintStream object to print stderr to
     * @param errMessage Error message to be added to error output (e.g. which cluster caused the error)
     * @param env A map of environmental variables to introduce. If null, no additional environmental variables will be introduced.
     * @throws IOException
     * @throws InterruptedException 
     */
    public static void runProcess(String programCommand, List<String> arguments, PrintStream out, PrintStream errOut, String errMessage, Map<String, String> env) throws IOException, InterruptedException{
        List<String> command = new ArrayList<>();
        command.add(programCommand);
        if (arguments != null){
            command.addAll(arguments);
        }
        try{
            String line;
            ProcessBuilder builder = new ProcessBuilder(command);
            if (env != null){
                for (Map.Entry<String, String> entry : env.entrySet()){
                    builder.environment().put(entry.getKey(), entry.getValue());
                }
            }
            final Process process = builder.start();
            BufferedReader brCleanUp = new BufferedReader(new InputStreamReader(process.getInputStream()));
            while ((line = brCleanUp.readLine()) != null) {
                if (out != null) {
                    out.println("[" + programCommand + " stdout]: " + line);
                }
            }  
            brCleanUp.close();
            
            brCleanUp = new BufferedReader(new InputStreamReader(process.getErrorStream()));
            while ((line = brCleanUp.readLine()) != null) {
                if (errOut != null){
                    errOut.println(errMessage + "[" + programCommand + " stderr]: " + line);
                }
            }
            process.waitFor();           
            brCleanUp.close();
        } catch (IOException e){
            throw new IOException(e);
        }     
    }
}
