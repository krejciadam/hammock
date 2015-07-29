HAMMOCK - peptide clustering

Runnable version, manual and more info at: http://www.recamo.cz/en/software/hammock-cluster-peptides/

Galaxy version at: https://toolshed.g2.bx.psu.edu/view/hammock/hammock


BUILDING FROM SOURCE

Hammock is a java application, so there should be no need to build it from source, unlike external
packages used (Clustal Omega, Hmmer3, HHsuite), that may need building on local machine. Refer to Hammock
manual for more info. 

The Hammock.jar file, which is part of the runnable version, should run on any system with Java 
(version at least 1.7.0) and there should be no need to build it from source. If for some reason
there is a need to build Hammock java code from source, it can be done this way: 

- suppose we are in a directory containing folder "src", which contains a subfolder called "hammock"
containing all the .java source files. On a Unix-like machine, run these commands: 

javac src/hammock/*.java
cd src
echo "Main-Class: hammock.Hammock" > manifest.mf
jar cmvf manifest.mf ../Hammock.jar hammock/*.class

This generates a file called "Hammock.jar" in the original folder. To run Hammock using the code 
just compiled, download runnable version from (http://www.recamo.cz/en/software/hammock-cluster-peptides/)
or a release from the github page (https://github.com/hammock-dev/hammock/releases), extract files from the 
archive and in folder "dir", replace the original "Hammock.jar" file with the newly built file.





