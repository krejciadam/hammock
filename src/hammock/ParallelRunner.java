/*
 * Class performs parallel execution of Callables.
 */

package hammock;

import java.util.Collection;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;

/**
 *
 * @author Adam Krejci
 */
public class ParallelRunner {
    
    /** Executes a Collection of Callables in parallel and waits
     * for their completion.
     * 
     * @param callables Callables to be executed
     * @param threadPool ExecutorService to be used
     * @throws InterruptedException
     * @throws ExecutionException 
     */
    public static void ExecuteInParallel(Collection<Callable<Void>> callables, ExecutorService threadPool) throws InterruptedException, ExecutionException, Exception{
        CompletionService<Void> resultPool = new ExecutorCompletionService<>(threadPool);
        for (Callable<Void> callable : callables){
            resultPool.submit(callable);
        }
        for (int i = 0; i < callables.size(); i++){
            resultPool.take().get();
        }
    }
    
}
