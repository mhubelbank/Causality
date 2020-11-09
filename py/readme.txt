# To run various commands

cd <analytics_active>/py

# Run the system test harness with an optional test name
# Note: Test names are defined at the top of the file IC_TEST.py
python IC_TEST.py [<testname>]

# Plot each variable along with its estimated density (PDF)
# Note: Testfile can be .py or .csv.  It is the full path e.g. ..\tests\synthTest1.py
python synthDataPlot.py <testfile>

# Generate synthetic data from a SEM
python synthDataGen.py <testfile> [<datapoints>]

# Discover parameterised SEM from given data
python modelTest.py <testfile> [<datapoints>]

# Assess the power of various Causal Direction Tests against a pair of variables from a designated data file (.py or .csv).
# Run repeatedly <runs> times to assess the test's power against different generations of the data 
python directionTest.py <testfile> <var1> <var2> <datapoints> <runs>

# Latest run = python attractorPlot2.py ../Tests/attractor.py s1 s3 2000
