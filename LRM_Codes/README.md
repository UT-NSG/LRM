This folder contains the latest version of the codes that generate the results of the paper. 
Using the codes provided in this folder is recommended. 
Main file that you would need to run is "IterOver_MainRunningFile.R". All the parameters and settings are set inside this file. 
To run the code:
1. Make sure the working directory is set currectly to the directory the codes are in. (See line 18).
2. Make sure you have installed the required libraries, especially LinkPrediction. (See https://github.com/galanisl/LinkPrediction/blob/master/README.md for more information).
3. Make sure you have downloaded all required functions and have put them in their corresponding directories.
4. Make sure you have all the required datasets.
5. Change the datasets you want to test in the file "data_spec.txt".

Different parts of the code have some comments that can help understanding the code. 
The following description also tries to explain how to use the main running file of the code. (The following description will be completed over time). 
Description of parameters and setting in the "IterOver_MainRunningFile.R":
1. NumberOfLeadingEigenValues_Vector is a vector that determines the different eigenvalues the methods will use.
For example, NumberOfLeadingEigenValues_Vector = c(1, 5, 10, 15, 20, 25, 30, NA)
NA means all. Obviously, the number of eigenvalues cannot exceed the number of nodes in the network. Also, if it is so close, the method could be numerically unstable. 

2. 
