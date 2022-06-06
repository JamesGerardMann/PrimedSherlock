# PrimedSherlock

A python-based command-line script to augment CRISPR-Cas12 and CRISPR-Cas13 Assay generation using RPA Primer Sets. 


### Installation

```
gh repo clone JamesGerardMann/PrimedSherlock
```

### Parameters
```
MismatchThreshold
```
This sets the tolerance for generated crRNA target to mismatches present in the strains utilized by the tool. It should be somewhere between 0-3. Most of the examples utilized 3, where as the ZIKA example utilized 5. If no crRNA sets are generated you may want to increase this value.


```
Use_PrimedRPA
```
This sets the tool into a mode where it seeks a generated PrimedRPA Output_Sets.csv file. PrimedRPA is the best tool to generate primer sets for our tool. If you set this to No, you will need to provide hand generated primer pairs and ensure they are accurate. 


```
Threadcount
```
This sets the amount of CPU threads to use for Cas-Offinder instances. More equals faster data processing.



```
Amplicon_Max_Size
```
Sets limits on the max size of amplicon that the tool will use. This is important as larger amplicons equals slower amplification by RPA polymerase. 



```
Amplicon_Min_Size
```
Sets lower limits on the size of amplicon that the tool will use. This is important as smaller amplicons equals faster amplification by RPA polymerase. 

```
Max_Background_Identity
```
This is a background identity limit introduced by PrimedRPA. Lower is better, sets a off-target amplification threshold. Hypothetically shouldn't matter with how our tool generates crRNA instead of probes. 


```
Off_Target_Mismatch_Score
```
This sets lower limits on how many mismatches are required for the crRNA to be deemed safe from recognition of off-target sequences. 


```
Contingency_crRNA_Utilization
```
This is very experimental and should not be used. It recognizes which strains have low affinity to be detected by crRNA targets and generates potential seperate pairs for highly conserved, highly mismatched strains. 


```
Working_Directory, str(master_dir))
Primed_RPA_Output, str(master_dir))
Input.fna_Location, str(master_dir))
Background_Blast_File, str(master_dir))
```
This is used for developer troubleshooting. This is automatically set when the code first run. 


```
Run_Cas-off
```
This is used for developer troubleshooting. Turning this off will stop the code prior to the Cas-Offinder.exe step. Keep this on. 


### Key Output Files

PrimedSherlock runs generate the following output file:

```
[RunID]_FinalOutput.csv
```


### 3rd-Party Software

PrimerSherlock incorporates segments of ...


PrimedRPA - https://github.com/MatthewHiggins2017/bioconda-PrimedRPA


PrimedSherlock includes an excutable from...


Cas-Offinder http://www.rgenome.net/cas-offinder/

### Contact

If you encounter any bugs please contact me directly at **James_Mann1[at]alumni.baylor.edu**
