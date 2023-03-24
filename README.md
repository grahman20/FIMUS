# FIMUS
FIMUS imputes numerical and categorical missing values by using a data set’s existing patterns including co-appearances of attribute values, correlations among the attributes and similarity of values belonging to an attribute.

# Reference

Rahman, M. G. and Islam, M. Z. (2014): FIMUS: A Framework for Imputing Missing Values Using Co-Appearance, Correlation and Similarity Analysis, Knowledge-Based Systems, 56, 311 - 327, ISSN 0950-7051, DOI information: http://dx.doi.org/10.1016/j.knosys.2013.12.005. 

## BibTeX
```
@article{rahman2014fimus,
  title={FIMUS: a framework for imputing missing values using co-appearance, correlation and similarity analysis},
  author={Rahman, Md Geaur and Islam, Md Zahidul},
  journal={Knowledge-Based Systems},
  volume={56},
  pages={311--327},
  year={2014},
  publisher={Elsevier}
}
```

@author Md Geaur Rahman <https://csusap.csu.edu.au/~grahman/>
  
# Two folders:
 
 1. FIMUS_project (NetBeans project)
 2. SampleData 
 
 FIMUS is developed based on Java programming language (jdk1.8.0_211) using NetBeans IDE (8.0.2). 
 
# How to run:
 
	1. Open project in NetBeans
	2. Run the project

# Sample input and output:
run:
Please enter the name of the file containing the 2 line attribute information.(example: c:\data\attrinfo.txt)

C:\SampleData\attrinfo.txt

Please enter the name of the data file having missing values: (example: c:\data\data.txt)

C:\SampleData\data.txt

Please enter the name of the output file: (example: c:\data\out.txt)

C:\SampleData\output.txt


Imputation by FIMUS is done. The completed data set is written to: 
C:\SampleData\output.txt