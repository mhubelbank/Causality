# Causality Framework

The HPCC Systems Causality Framework is currently under development.  This repository provides a Python module that implements all of the relevant algorithms for causal model representation, validation, causal inference, and ultimately counterfactual inference.  It incudes a synthetic data generation capability, used for testing the algorithms, and for probing their power against various scenarios.  This repository can be used as a stand-alone framework.  It will also be used to provide local processing
for a parallelized causality framework on HPCC Systems supercomputing clusters.
For details on Causality, Causal Inference, and the various causality algorithms, please see the References section below.

## Installation
### Prerequisites
- Python3
- numpy
- matplotlib
### Procedure
- Clone repository
- Add your repository to pythonpath or work within the repository.

## Causal Methods

_cGraph.py_ (Causal Graph) is the heart of the system.  It accepts a multivariate dataset as well as a Causal Model (aka Path Diagram) that is thought to correspond to the causal process that produced the data.  The Causal Model is a Directed Acyclic Graph that represents the best hypothesis of the causal process that underlies the generation of the provided dataset.  The Causal Model may contain "unobserved" variables as well as the "observed" variables with their corresponding datasets.  It includes the following functional areas:
- TestModel() analyzes the data and highlights areas where the actual data deviates from the Causal Model's Hypotheses.  It enumerates the independence implications of the Causal Model, and verifies each of these expected independencies / dependencies using various conditional independence testing methods.
- Intervene() implements the "Do Calculus", for causal inference, emulating the results of a controlled experiment.  The Do Calculus attempts to screen-off non-causal correlation and identify the effect of an intervention (i.e. setting a variable to a fixed value) on other variables. If the model does not provide sufficient observed variables to screen-off the non-causal effects, than the result is returned as "not identifiable".  Intervention is the method that implements.
- Causal Metrics -- Various Causal Metrics can be requested regarding a pair of variables including: Average Causal Effect (of X on Y), Controlled Direct Effect, and Indirect Effect.
- Counterfactual Inference (Future) -- Provides a method for counterfactual inference such as "Would John have recovered from his illness if he had been given a certain treatment, given that he was not given the treatment and did not survive".
- Model Discovery (Future) -- Attempt to infer a Causal Model from the data, by deep analysis of the provided data.

## Conditional Probability

_Probability/Prob.py_ provides a powerful multivariate conditional probability analysis system.  It equally supports both discrete (including categorical) and continous data, or any mixture of the above.  Given a set of data, Prob allows a wide range of queries including:
- P(X=x) -- The numeric probability that X takes on the value x.
- P(X=x, Y=y) -- Joint Probability -- The probability that X=x and Y=y simultaneously.
- P(Y=y | X=x) -- The numerical probability that Y takes on the value y, given that X is observed to have the value x.
- P(Y | X=x, ... , Z=z) -- The conditional probability distribution of Y, given the values of multiple variables.
- P(Y=y | X=x, ... , Z=z) -- The numerical probability that Y takes on the value y, given the values of multiple variables.
- P(Y=y | X) -- The numerical probability that Y=y when conditionalized on X.  This is equivalent to the sum of P(Y=y | X=x) * P(X=x), for every value x of X.
- P(X) -- The probability distribution of any variable X.  Presented as a PDF (Probability Distribution Function), which is actually a discretized version of the probability distribution.
- P(Y | X=x) -- Conditional Probability -- The probability distribution of Y given that X is observed to have the value x.
- P(Y | X=x, ... , Z=z) -- The probability distribution of Y, given the values of multiple variables.
- P(Y | X) -- The probability distribution of Y when conditionalized on X.  This is equivalent to the sum of P(Y | X=x) * P(X=x), for every value x of X.
- P(Y | X=x, Z) -- The probability distribution of Y given that X=x and conditionalized on Z.  There may be additional givens, and multiple variables may be conditionalized on.
- P(Y=y | X=x, Z) -- The numerical probability that Y=y,  given that X=x and conditionalized on Z. There may be additional givens, and multiple variables may be conditionalized on.
- P(X=x, Y=y, ... ,Z=z) -- The numerical joint probability of multiple events occurring simultaneously.
- Dependence(X, Y) -- The numerical dependency measured between X and Y.
- Dependence(X, Y | Z=z, W=w) -- The numerical dependency measured between X and Y given one or more observed covariate values.
- Dependence(X, Y | Z=z, W) -- The numerical dependency measured between X and Y, given that Z=z and conditionalized on W.  There may be multiple conditional values or variables to conditionalize on.

This module also

Probability/pdf.py Provides a Probability Distribution Function object, which is a finely discretized probability distribution.  This is returned by Prob (above) whenever a probability distribution is requested.  The following functions are available for a PDF:
- pdf.P(v) -- The probability of a value "v" occurring.
- pdf.E() -- The Expected Value of the distribution (i.e. mean).  This is the first "standard moment" of the distribution.
- pdf.stDev() -- The Standard Deviation (sigma) of the distribution. This is the second standard moment of the distribution.
- pdf.skew() -- The skew (third standard moment) of the distribution.
- pdf.kurtosis() -- The excess kurtosis (fourth standard moment) of the distribution.
- pdf.compare(otherPdf) -- Compare two PDFs, and assess their level of similarity in terms of the four moments.

## Synthetic Data Generation

The "models" folder contains various test models including a set of standard test models M0 - M8.  These are synthetic data models, each of which contains a Structural Equation Model (SEM) that is used to generate a dataset with known causality, as well as a Causal Model that represents an hypothesis as to the data generation process.  Generally these two aspects are set up to be consistent, but can be made inconsistent to e.g., see the result of the model validation process.  The SEM can be linear or non-linear, and can include any mathematical functions, or noise distributions (e.g. Normal, Logistic, Exponential, Beta, LogNormal).  The Causal Model may contain variables not described by the SEM, in which case, those variables are treated as "unobserved".  The command line utility synth/synthDataGen.py (see below) is used to process the SEM and generate datasets of any size.  These are produced as CSV datasets.  The synth/getData.py module is used to read these CSV datasets.

### Using synthDataGen:

    python3 synth/synthDataGen.py <model file> <numRecords>

    Example:
        python3 synth/synthDataGen.py models/M3.py 100000
        Genrates 100,000 multivariate recordds as described by M3.py

### Assessing a synthetic model:

    The command line program cTest will utilize Graph.py and a synthetic model with generated data to determine the extent to which the defined Causal Model is consistent with the generated data.

    python3 cTest.py <model file> [<numRecords>]

    Example:
        python3 cTest.py models/M3.py
        python3 cTest.py models/M3.py [10000]


    If <numRecords> is not provided then all generated records will be used.
    If <numRecords> is provided, it will take a sample of <numRecords> from the generated data.  
    If <numRecords> is greater than or equal to the number of generated records, then all records will be used.


## Other command-line programs:

interventionTest.py provides a test and demonstration of the intervention logic, and causal metrics.  It is hard coded to test between variables 'A' and 'C', but can be modified locally to test any desired variables.  It demonstrates use of cGraph and  getData.

    python3 interventionTest.py <model file>

    Example:
        python3 interventionTest.py models/M3.py

probTest.py is a test for the probability system, and demonstrates the use of the probability system using synthetic data.  It is hardcoded for particular models but can be modified locally for testing any model.

    python probTest.py <model file>

    Example:
        python3 probTest.py models/M1.py

## References:

[1] Mackenzie, D. and Pearl, J. (2018) _The Book of Why._ Basic Books.

[2] Pearl, J. (2000, 2009) _Causality: Models, Reasoning_, and Inference.  Cambridge University Press

[3] Pearl, J., Glymour, M., Jewell, N. (2016) _Causal Inference In Statistics._ John Wiley and Sons.

[4] Morgan, S.L. and Winship, C. (2014) _Counterefactuals and Causal Inference: Methods and Principles for Social Research, Analytical Methods for Social Research._ 2nd edn.  Cambridge University Press.

[5] Shimizu, S. and Hyvarinen, A. (2006) _A linear non-Gaussian aclyclic model for causal discovery._ Journal of Machine Learning Research, 7:2003-2030.

[6] Spirites, P. and Glymour, C. (1991) _An algorithm for fast recovery of sparse causal graphs._ Computer Review, 9:67-72.


