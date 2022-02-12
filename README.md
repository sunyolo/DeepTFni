# DeepTFni
A VGAE-based model to infer transcription factor regulatory network (TRN).
TRN is represented as an undirected graph, where V nodes represent TFs and 
E edges represent their interactions. The input of DeepTFni is solely a scATAC-seq
count matrix. The output of DeepTFni is the imputed TRN.

## Workflow
TRN is represented as an undirected graph, where V nodes represent TFs and 
E edges represent their interactions. The input of DeepTFni is solely a scATAC-seq
count matrix. The output of DeepTFni is the imputed TRN. Taking TRN inference as a link prediction task, DeepTFni workflow is organized as follows:
![image](https://github.com/sunyolo/DeepTFni/blob/main/data_resource/DeepTFni%20workflow.pdf)

## Dependency
DeepTFni is mostly written in Python 3.6 and some preprocessing steps are done in Perl 5. It can be run on a single desktop using Linux platform. To run DeepTFni, some prerequisites need to be installed. A detailed dependency list is:
* pytorch 1.7.1
* numpy 1.8.11
* pandas 1.1.3
* matplotlib 3.3.2
* MAESTRO 1.5.0
* BioPerl 1.5.7
* R 4.0.3
* fimo

## Running
Use DeepTFni through command:
```
perl run_DeepTFni_csv.pl Path_to_Input reference #Bash
``` 
For example:
```
perl run_DeepTFni_csv.pl ./Input_data/Input.csv hg19 #Bash
``` 
It will generate a new script based on run_template_csv.sh named like ‘run_your_input.sh’. 
The input should be organized like ./demo/scATAC_demo.csv, that is, each row represents a peak and each column represents a cell. If DeepTFni goes smoothly, it will finally return the imputed TRN represented by a symmetric adjacency matrix stored at ./train_info/result_k_10_summary. 
