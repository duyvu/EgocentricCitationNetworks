# EgocentricCitationNetworks

This C++ projects implements egocentric event models for citation networks.

The corresponding research paper was published in Proceedings of the 28 th International Conference
on Machine Learning, Bellevue, WA, USA, 2011. A copy of the paper is stored in the directory **papers**.


## Compiling the code:

1. Unzip the file

2. Go to the unzipped directory

3. Download boost library from http://www.boost.org/users/download/ and unzip the file. No compilation for boost library is needed.

4. Specify directories in the file "makefile":

* HBOOST : path to where you unzip the boost library

* HRLIB : path to R header files, e.g. /usr/share/R/include

* HMYDIR : path to the source folder

* RLIB = path to the directory containing compiled libraries of R, e.g. /usr/lib/R/lib

5. Run "make"

## Running the code:

We assume that the data set arXiv-TH is placed in the directory "data/CitationHepTh" (I include them in the zip file already). There are 3 files for nodes, edges, and LDA covariates. The paths to these files are hard code in icmlArXivHepTh.h but can be easily modified.

To run the code just type

LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH:lib
export LD_LIBRARY_PATH
./egonet modelType startTime endTime signature

The first two lines are needed only in case you have not set the lib path:

1. modelType can be 0 to 5. For example 0 is for preferential attachment model only, 4 is for LDA only, 3 is a combination of LDA and other network statistics (the best model in the paper). You can check icmlArXivHepTh.h to see other options.

2. startTime and endTime are used to specify the time interval where events are used to fit the Cox model. Events before the startTime are still used to construct the network statistics. Events after endTime are just ignored.

3. You can put any text into signature. This text will add to the end of any output file so that you can discriminate outputs of different experiments.

One example is 

./egonet 0 11464 11585 1000

It will fit a PA model using events from the day 11464 (since 1970) to 11585. The output file have a suffix 1000 and is placed in the directory output. You can look at the screen to see Newton-Raphson iterations and the final estimated coefficients. The output file contains information about the model type, start and end times, the likelihood at MLE, coefficients, and covariance matrix.
