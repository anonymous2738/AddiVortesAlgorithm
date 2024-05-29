AddiVortes
===========

An R (Bayesian) Additive Voronoi Tessellation implementation (AddiVortes)
Software for Supervised Statistical Learning

Copyright (C) 2024
Adam Stone  
Department of Mathematical sciences, Durham University
& 
John Paul Gosling
Department of Mathematical sciences, Durham University
 
Setup Instructions
------------------

To install the AddiVortes function in R, you can run the following code: 

```r
source_url("https://raw.githubusercontent.com/anonymous2738/AddiVortesAlgorithm/main/AddiVortesMainCode.R")
# Source the raw script from GitHub
source_url("https://raw.githubusercontent.com/anonymous2738/AddiVortesAlgorithm/main/Datasets.R")

source_url("https://raw.githubusercontent.com/anonymous2738/AddiVortesAlgorithm/main/CodeForFigures.R")
```


### Install Java JDK (not the JRE)

Download the latest [Java JDK](https://jdk.java.net/) and install it properly. (Java 7 or above is required for v1.2.x and Java 8 or above is required for v1.3 and above). bartMachine requires rJava which requires the JDK; you cannot just have a JRE!

### Install rJava

Use `install.packages("rJava")` within R. If you experience errors, make sure your `JAVA_HOME` system variable is set to the root of your java installation (on a windows machine that would look something like `C:\Program Files\Java\jdk-13.0.2`). Also try running `R CMD javareconf` from the command line. On ubuntu, you should run `sudo apt-get install r-cran-rjava` to install from the command prompt. If you still have errors, you are not alone! rJava is tough to install and idiosyncratic across different platforms. Google your error message! The majority of issues have been resolved on Q&A forums.

### Install bartMachine via CRAN

Use `install.packages("bartMachine")` within R.

### Install bartMachine via compilation from source

Due to CRAN limitations, we cannot release bartMachine for Java >7. Thus it is recommended to install bartMachine from source. This is recommended as you will get the benefits of Java 8-14 as well as the latest release of bartMachine if it's not on CRAN yet (see [changelog](https://github.com/kapelner/bartMachine/blob/master/bartMachine/CHANGELOG).

1. Make sure you have [git](http://git-scm.com/downloads "Download git for all operating systems") 
properly installed.

2. Run `git clone https://github.com/kapelner/bartMachine.git` from your command line and navigate into the cloned project directory via `cd bartMachine`.

3. Make sure you have a [Java JDK](https://jdk.java.net/14/) installed properly. Then make sure the bin directory is an element in the PATH variable (on a windows machine it would look something like `C:\Program Files\Java\jdk-13.0.2\bin`). We also recommend making a system variable `JAVA_HOME` pointing to the directory (save \bin).

3. Make sure you have [apache ant](http://ant.apache.org/bindownload.cgi "Download apache ant for all operating systems") installed properly. 
Make sure you add the bin directory for ant to your system PATH variable (on a windows machine it would be something like `C:\Program Files (x86)\apache-ant-1.10.8\bin`). We also recommend making a system variable `ANT_HOME` pointing to the directory (save \bin).

4. Compile the JAVA source code into a JAR using `ant`. You should see a compilation record and then `BUILD SUCCESSFUL` and a total time.

5. Now you can install the package into R using `R CMD INSTALL bartMachine`. On Windows systems, this may fail because it expects multiple architectures. This can be corrected by running `R CMD INSTALL --no-multiarch bartMachine` (I haven't seen this issue in years though). This may also fail if you don't have the required packages installed (run `install.packages("bartMachineJARs")` and `install.packages("missForest")`). Upon successful installation, the last line of the output should read `DONE (bartMachine)`. In R, you can now run `library(bartMachine)` and start using the package normally.


#### Limiting CPU usage
(At least under GNU/Linux) even if you set `set_bart_machine_num_cores(1)`, CPU usage per process can be much larger than 100% (reaching at times 200% or 300%). This can lead to CPU overloading, especially if you run multiple bartMachines in parallel (for example, if you use the [SuperLearner](https://cran.r-project.org/web/packages/SuperLearner/) package and use parallelization). This seems to be a consequence of the garbage collector. One way to avoid this problem is to issue `Sys.setenv(JAVA_TOOL_OPTIONS = "-XX:ParallelGCThreads=1")` *before* invoking `library(bartMachine)`. (If you use a cluster, for example a SNOW cluster, you will want to do this in the slaves too, for example `clusterEvalQ(the_name_of_your_cluster, {Sys.setenv(JAVA_TOOL_OPTIONS = "-XX:ParallelGCThreads=1")})`).

Acknowledgements
------------------

We thank Ed George, Abba Krieger, Shene Jensen and Richard Berk for helpful discussions. We thank Matt Olson for pointing out an important memory issue. We thank [JProfiler](http://www.ej-technologies.com/products/jprofiler/overview.html) for profiling the code which allowed us to create a lean implementation.
