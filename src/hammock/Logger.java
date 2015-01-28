/*
 * Writes information on stderr and to run.log file
 */
package hammock;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.Calendar;

/**
 *
 * @author Adam Krejci
 */
public class Logger {

    private final String filePath;
    private final SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss.SSS");
    private Calendar calendar;
    private final boolean dummyLogger;

    public Logger(String filePath, boolean dummyLogger) {
        this.filePath = filePath;
        this.dummyLogger = dummyLogger;
    }

    /**
     * prints line to output file only
     *
     * @param line
     */
    public void logWithoutTime(String line) {
        if (dummyLogger) {
            return;
        }
        if (filePath != null) {
            try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(filePath, true)))) {
                writer.println(line);
            } catch (IOException e) {
                System.err.println("Warning: Failed to log follwoing message. Run will continue, message will not be appended run.log");
            }
        }
    }

    /**
     * Adds time stamp to line and prints it to output file only
     *
     * @param line
     */
    public void logWithTime(String line) {
        if (dummyLogger) {
            return;
        }
        calendar = Calendar.getInstance();
        String newLine = dateFormat.format(calendar.getTime()) + ":\t" + line;
        logWithoutTime(newLine);
    }

    /**
     * Adds time stamp to line and prints it to output file and stderr.
     *
     * @param line
     */
    public void logAndStderr(String line) {
        if (dummyLogger) {
            return;
        }
        logWithTime(line);
        System.err.println(line);
    }
}
