# gp_general program
***
## Description
***

The gp_general program is a Fortran90 implementation of multidimensional
Gaussian Process that can learn from value and/or derivative observations
with or without sparsification and predict function values, first and second 
derivatives and corresponding variances. Covariance can be calculated based 
on the squared exponential kernel that is hard coded and massively parallelized 
or arbitrary user defined kernel. The implementation uses OpenMP multithreading.

***
## Authors
***

L. Mones, N. Bernstein and G. Csanyi

***
## Get the code
***

For a more convenient code development/maintenance the program package is divided into two parts ([prerequisites](https://github.com/molet/prerequisites) and [project_gp_general](https://github.com/molet/project_gp_general)). Both are available on GitHub. To get all components just simply run the following script:

	./Pull_code

***
## Installation
***

 To install the software use cmake in the main directory and cross your fingers:

	cmake .

If you are lucky the build process is successful and you have a `Makefile`. Then
simply type:

    make

The executable will be installed into the bin directory (`gp_general`).
Both intel and gfortran compilers should work. If multiple compilers available
use the standard way in cmake to specify the one you want (i.e. supposing you are using
bash: `export FC=gfortran)`.
To use preinstalled BLAS/LAPACK libraries you have the following options:
* Intel MKL library: if environment variable MKLROOT is specified then it is
  automatically detected and used (enforce this option use `-DBLAS_VENDOR=MKL`)
* Apple's Accelerate framework for MacOSX: use `-DBLAS_VENDOR=ACC`
* Specifying path to BLAS/LAPACK libraries: use `-DBLAS_LAPACK_LIBRARIES=`/path/to/lib
  (enforce this option by using `-DBLAS_VENDOR=LIB`)
If no preinstalled libraries of the above specifications are found then internal BLAS
and LAPACK routines provided by the package are compiled and applied. This routines 
are not optimised so to achieve good performance for double precision calculations
internal BLAS/LAPACK is not advised. 
Quad precision is also available, for this type:

	cmake -DPREC=QUAD

This will use the internal linear algebra subroutines.
By default subprojects are built as shared libraries. Using them as static add
`-DLIB_BUILD_TYPE=STATIC`. 
Debug mode can be activated by the option `-DCMAKE_BUILD_TYPE=DEBUG`.

***
## Usage
***

 The list of available topics can be obtained by simply calling the executable without
any argument or with the standard help flags:

	gp_general
        gp_general -h[elp]
        gp_general --help

Calling a given topic will list the corresponding options.
Options must be defined as arguments after `gp_general`. The arguments are case-insensitive
and the order is arbitrary.
Below there is a short explanation of the possible options.

***
### Specifying the type of input
***

* The program can learn from values (`-InVal`), derivatives (`-InDer`) or both (`-InVal -InDer`).
* The data input file's name can be specified by `-TrainFile` *trainfile*.
By default the program looks for a file called `training.txt`. The format of the data file is
discussed later (see **File formats** section).
* In the case of sparsification, previously taught data that are saved in one ore more gpfiles can be also used by specifying 
`-InGP`. The files should be specified as `-GPInFile` *gpinfile1* `-GPInFile` *gpinfile2* etc. 
If no input data file is specified then the gpfiles are simply merged together, otherwise the new teaching result is also added to 
the previous one(s).

***
### Specifying input variables
***

* The number of variables has to be defined by `-NDim` *ndim*. This will allocate several arrays
that are used even in the input parsing process.
* `-Periodicity` *per1* *per2* etc. defines the possible periodicities of the variables 
(defining 0 periodicity means that the variable is not periodic).
* The available interval of interest is either calculated automatically
from the actual training points using the minimum and maximum values for each component or can
be specified by `-Range` *min1* *max1* *min2* *max2* etc. Where a grid is required by default this range
is used unless the grid is explicitely defined (see `-PredRange`, `-SparseRange` and `-MinRange`).
* By default all lines of the input datafile is used but this can be controlled by the 
`-LStart` *lstart* (first line to read), `-LStop` *lstop* (last line to read) and `-LInc` *linc* (frequency of lines to read) keywords.
* The first column of the datafile can be omitted by applying `-RowNum`.
* There are two types of kernel that can be specified by the keyword `-Kernel`: squared exponential
(`se`) and user defined (`user`).
* In the case of `se` kernel the function deviation (`-Delta` *delta*),
the theta values (`-Theta` *theta1* *theta2* etc.) must be specified. Function values can be shifted
by `-AveVal` *aveval* before the training.
* User defined kernel must be specified by `-FKernel` "...". Please use quotation marks (single
or double) if you apply white space characters in this definition. The function synthax must 
be standard fortran synthax. First and second variables have the form of x1, x2, etc and 
y1, y2, etc. where the numbers denote the corresponding components (up to *ndim*). Using 
`-Prefix1` *prefix1* and `-Prefix2` *prefix2* the x and y prefixes can be modified (no white space 
character is allowed for prefix definitions).
* Value and derivative deviations can be provided either in the input datafile (e.g. if they are 
known and different) or uniform values can be specified by `-SigmaVal` *sigmaval* and `-SigmaDer` *sigmader1*
*sigmader2* etc.
* `-Jitter` *jitter* adds a small number to the diagonal elements of the covariance matrices to get
a better conditioned matrix.
* For the calculation of the inverse of the covariance matrix several factorization methods are
available. By default the Cholesky factorization is used but this can be modified to QR (or
Bunch-Kaufman) using `-Factorization` `qr` (or `bk`).
* By default no sparsification is applied. The sparsification can be activated by using the
-Sparsification with a corresponding keyword. There are three different ways currently implemented to pick the
sparse points: 
  * `clustering` uses a cluster analysis based on the provided data points. The clustering method
   can be specified by the `-Clustering` keyword that can be `every` (sparse points equidistantly
   selected from the data file), `kmeans` (using kmeans clustering method) and `kmeans++` (using 
   kmeans++/kmeans clustering method). The number of sparse points is defined by `-SparsePoints`
   *sparsepoints*. For kmeans and kmeans++/kmeans methods the option `-SparseFreq` *sparsefreq* can be 
   used to apply the methods to only a subset of data points using every *sparsefreq*-th point. 
   `-KMeansThreshold` can be used to accelerate kmeans method by specifying thresholds for each 
   component.
  * `grid` uses an equidistant grid based on the interval calculated automatically or specified 
   by `-Range`. If different range is needed then use `-SparseRange` *min1* *max1* *min2* *max2* etc. 
   The grid must be specified using the `-SparseGrid` *grid1* *grid2* etc. keyword.
  * `file` uses an external file that includes the sparse points positions. The name of the file
   can be given by `-SparsePosFile` *sparseposfile* and if the first column is needed to be skipped the
   keyword `-SparseRowNum` must be used.
* Based on the article "A Unifying View of Sparse Approximate Gaussian Process Regression" by J.
Quinonero-Candela and C. E. Rasmussen (Journal of Machine Learning Research 6 (2005) 1939) the
following sparsification techniques can be used: (`-SparseMethod`) `dic`, `dtc` (default) and `fitc`.
The type of sparse points can be selected using the `-SparseType` option and it can be value 
(`val`, default), derivative (`der`) or both (`all`).
* Random seed for random number generator can be specified by using `-Seed` *seed*.
* `-Verbose` activates a verbose mode.

***
### Specifying the type of output
***

* The program can predict values (`-OutVal`), their deviations (`-OutValDev`), derivatives
(`-OutDer`) and their deviations (`-OutDerDev`), and in the case of squared exponential kernel
even integrals (`-OutInt`) and their deviations (`-OutIntDev`).
* When integrals are prediceted integration intervals must be specified using `-Interval` *int_min1* *int_max1* *int_min2*
*int_max2* etc. Variables with the same *int_min* and *int_max* are not integrated out.
* Prediction points can be defined by either using `-PredGrid` *grid1* *grid2* etc. (which applies an
equidistant grid) or providing the points explicitely: `-PredPosFile` *predposfile*. Skipping the first
column of this file `-PredRowNum` has to be added.
* By default, results of prediction are written to `pred.txt` file. The name of this file can be
modified by using `-PredFile` *predfile*.
* When values are predicted, the final result can be shifted by `-Shift` *shift* to move the minimal value to 0.
* `-Gnuplot` keyword inserts empty lines between blocks if grid based prediction was requested. This is very useful for 3D visualization using gnuplot.
* For later use a special gpfile can be printed out (`-OutGP`). The name of this file
can be specified by `-GPOutFile` *gpoutfile*, by default it is `gp.txt`. Using gpfile is 
very useful when sparsification is applied for large data set by splitting the data file
to smaller portions, processing them separately and then merging them together.
* There are three different ways to store data in gpfile using the keyword `-GPMode`. The default is `full`
that prints out all information. This is useful for conventional GP but could be very
large when sparsification is used. For such situations the more economical `partial` 
mode is suggested. It also contains all the informations required for prediction but
the inverse of the covariance matrix is not calculated and stored in the gpfile. Finally,
`compact` mode results in a realtively small gpfile that includes only the information
for doing prediction but no more training data can be added.
* Hyper parameter optimization is available for squared exponential kernel (`-OutOpt`). `-OptDelta`, 
`-OptSigmaVal`, `-OptSigmaDer` and `-OptTheta` keywords specify the relative stepsizes for each parameter
and using 0 value does not optimize the corresponding parameter. The result of the optimization
is written into `opt.txt` file that can be altered by `-OptFile` *optfile* keyword. Numerical 
differentiation is used to calculate gradients in the parameter space. `-NumDiff` defines whether
it is `one-sided` (default) or `two-sided`. Optimization related keywords are explained in 
section **Optimization parameters**.
* Searching for local minima is also possible on the reconstructed GP surface (`-OutMin`). Walkers
are launched from the prediction positions and minimizations are performed. Specific grid can be 
specified by `-MinGrid` *grid1* *grid2* etc. or external file including the starting positions can be
used (`-MinPosFile` *minposfile*). Specific range can be defined by `-MinRange` *min1* *max1* *min2* *max2* etc.
The found minima are printed out into the `min.txt` file, the name of this file can be altered by 
using `-MinFile` *minfile*.

***
### Optimization parameters
***

* Optimizer can be selected by `-Optimizer` *optimizer*. Currently only limited memory BFGS method
can be used with bound-constraint.
* `-NCorr` *ncorr* specifies number of corrections for lbfgs optimizer. 
* Maximum number of iteration can be given by `-MaxIter` *maxiter*.
* `-MaxDVal` *maxdval*, `-MaxGrad` *maxgrad*, `-MaxRMSGrad` *maxrmsgrad* correspond to the specification of thresholds of 
maximum change of the function being optimized, maximum of the absolute value of gradient 
components and maximum RMS value of the gradients. If at least one of these criteria is satisfied,
the optimization stops.

***
### File formats
***

Below the required formats are specified for the different files. Each entry represents a column in
the file and the order is strict. Brackets mean that the entry can be skipped by appropriate
specification.

#### Training file format (input):
***

    [ROW_NUMBER] POS(1) POS(2) ... POS(ndim) [VAL] [VAL_VAR] [DER(1) DER(2) ... DER(ndim)] [DER_VAR(1) DER_VAR(2) ... DERV_VAR(ndim)]

#### Position file format for prediction (input):
***

    [ROW_NUMBER] POS(1) POS(2) ... POS(ndim)

#### Position file format for sparsification (input):
***

    [ROW_NUMBER] POS(1) POS(2) ... POS(ndim)

#### Position file format for minimization (input):
***

    [ROW_NUMBER] POS(1) POS(2) ... POS(ndim)

#### Prediction file format (output):
***

    POS(1) POS(2) ... POS(ndim) [VAL] [VAL_VAR] [DER(1) DER(2) ... DER(ndim)] [DER_VAR(1) DER_VAR(2) ... DERV_VAR(ndim)]

#### Hyper parameter optimization file format (output):
***

    [START_REL_FVAR=1] [START_REL_VALVAR=1] [START_REL_DERVAR(1)=1] [START_REL_DERVAR(2)=1] ... [START_REL_DERVAR(ndim)=1] \
    [START_REL_THETA(1)=1] [START_REL_THETA(2)=1] ... [START_REL_THETA(ndim)=1] -> \
    [OPT_REL_FVAR=1] [OPT_REL_VALVAR=1] [OPT_REL_DERVAR(1)=1] [OPT_REL_DERVAR(2)=1] ... [OPT_REL_DERVAR(ndim)=1] \
    [OPT_REL_THETA(1)=1] [OPT_REL_THETA(2)=1] ... [OPT_REL_THETA(ndim)=1] \
    NUM_OF_ITERATIONS -LOG_OF_LIKELIHOOD \
    [GRAD_REL_FVAR=1] [GRAD_REL_VALVAR=1] [GRAD_REL_DERVAR(1)=1] [GRAD_REL_DERVAR(2)=1] ... [GRAD_REL_DERVAR(ndim)=1]

#### Function minimization file format (output): 
***

    START_POS(1) START_POS(2) ... START_POS(ndim) -> OPT_POS(1) OPT_POS(2) ... OPT_POS(ndim) NUM_OF_ITERATIONS \
    FUNCTION_VALUE GRAD(1) GRAD(2) ... GRAD(ndim)
