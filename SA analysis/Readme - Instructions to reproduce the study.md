# S-aureus-MDR-antibiogram

This package reproduces all analysis found in "Analysis of multidrug resistance with a machine learning generated antibiogram for Staphylococcus aureus at one New York hospital between 2008 and 2018". 

Step 1.  Open CAZER R Project File 'S-aureus-MDR-antibiogram.Rproj', then getwd() to check that the path of this project file is the active directory.
Step 2.  Place data file "Linelist SENTRY site_136.xlsx" in the data folder.
Step 2.  In R, open the "master run.R" file and source/run it with echo.  "master run.R" calls all necessary scripts for the analysis.

VERY IMPORTANT NOTE: There is a possibility that the code will stop to prompt you to install Rtools and/or Java.  If this is the case, download and installation instructions are found here:
- Download and install RTools from here: https://cran.r-project.org/bin/windows/Rtools/
  Make sure to install folder to following location:  "C:\Program Files\R\Rtools"
- Download and install Java. Use the following link to download and install Java correctly: https://softfamous.com/jre7/
  After installing Java, if the code still prompts you for missing java, restart your computer or submit the following command into R before running the "master run.R" file:
  Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jre7')


Briefly, the analysis has the following steps:
1. "data prep.R": Prepare the data by interpreting MIC as susceptible or non-susceptible (also referred to as resistant)
2. "antibiogram and mdr profiles.R": Calculation of an antibiogram with the prevalence of resistance in the data overall and also stratified by year, MRSA vs MSSA, and infection type. Analysis of the multidrug resistance profiles (e.g. co-resistance to antimicrobial drugs and classes), including a simulation to demonstrate how the actual data varies from a null hypothesis of no association between resistances
3. "mining sets.R": Mining itemsets from the data stratified by year, MRSA vs MSSA, and infection type
4. "PValues.R": calculating the P-value of each itemset using a statitical simulation
5. "filtering and summarizing sets.R": selecting only itemsets with P-value <= 0.05 and collating them into one large dataframe
6. "bootstrap.R": Calculating the bootstrap 95%tile interval for itemset cLift
7. "AM circle plots.R" and "class circle plots.R": Visualization of itemsets in circle plots by decomposing itemsets into nodes (antimicrobials) and edges (associations)
8. "SST MRSA v MSSA.R":  a sub-analysis of MRSA and MSSA isolates from skin/soft tissue infections
9. "descriptive analysis.R": displaying descriptive statistics of specific itemsets used as examples in the manuscript
10. "methods example.R": verifying calculations used in Appendix tables to demonstrate association mining calculations

Data is available upon request and approval by the SENTRY Antimicrobial Surveillance Program. Aggregated data is available at https://sentry-mvp.jmilabs.com/app/sentry-public. License for use of this data will be determined by the SENTRY Antimicrobial Surveillance Program when data is released to investigators. The data file must be placed in the data folder.

Data Citation: Stephen Jenkins, JMI Laboratories. (2018). "Antimicrobial susceptibility testing SENTRY site 136 2008-2018". [Excel file "Linelist SENTRY site_136.xlsx"].

-------------------------------
 R.Version()
$platform
[1] "i386-w64-mingw32"

$arch
[1] "i386"

$os
[1] "mingw32"

$system
[1] "i386, mingw32"

$status
[1] ""

$major
[1] "3"

$minor
[1] "6.1"

$year
[1] "2019"

$month
[1] "07"

$day
[1] "05"

$`svn rev`
[1] "76782"

$language
[1] "R"

$version.string
[1] "R version 3.6.1 (2019-07-05)"

$nickname
[1] "Action of the Toes"


---

General information
---

Article Citation:

Title of the paper: Analysis of multidrug resistance with a machine learning generated antibiogram for Staphylococcus aureus at one New York hospital between 2008 and 2018
Journal:  

Author(s) and ORCID(s): 

- Casey L. Cazer - 0000-0002-2290-1868
- Lars F. Westblade -
- Matthew S. Simon - 
- Reed Magleby - 
- Mariana Castanheira - 0000-0002-7096-4310
- James G. Booth - 0000-0001-9572-6004
- Stephen G. Jenkins - 
- Yrjö T. Gröhn - 0000-0002-2721-9833

 Corresponding author: Casey L. Cazer
 Contact details:  email: clc248@cornell.edu 

Study Mnemonic: R2-2020-CAZER-1

---
Code Citation: 

Cazer, Casey L. (2020).   "Replication Package for Analysis of multidrug resistance with a machine learning generated antibiogram for Staphylococcus aureus at one New York hospital between 2008 and 2018" (version 1) [R files].  DOI:

Code Author ORCID:      Casey L. Cazer: 0000-0002-2290-1868
Code last updated: 	March 27, 2020	
Software and version:   RStudio Version 1.0.136; R 3.6.1
						Run on Windows 10, 64-bit

Code License: BSD-3-Clause

Copyright 2020 Casey L. Cazer
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer 
   in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from 
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

