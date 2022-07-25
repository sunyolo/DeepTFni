#! /bin/bash

sample_type="template"
echo "processing $sample_type"

reference="template"
echo "using reference $reference"

##### 	generate initial adjacency matrix	#####
command_1=`perl 1-ATAC_peaks.pl $sample_type $reference`
echo $command_1

#####	generate RP score using MAESTRO	#####
command_2=`MAESTRO scatac-genescore --format h5 --peakcount /Your_path/$sample_type/8_$sample_type\_peak_count.h5 --genedistance 10000 --species GRCh38 --model Enhanced --outprefix /Your_path/$sample_type/9_$sample_type`
#make sure to use Absolute Path for MAESTRO scatac-genescore function
echo $command_2

#####	    extract TF RP score as feature matrix	#####
command_3=`Rscript 7_for_TF_rp_score.R $sample_type`
echo $command_3

#####	        split data set for DeepTFni		#####
# require source ~/.bashrc
command_4=`python 8_data_preparation_for_DeepTFni.py -sample $sample_type`
echo $command_4

#####	               train DeepTFni	        	#####
command_5=`python 9_train_DeepTFni.py -sample $sample_type`
echo $command_5

#####	               organize output	#####
command_6=`python 10_weighted_2_binary.py -sample $sample_type`
echo $command_6
command_7=`python 11_integrated_output.py -sample $sample_type`
echo $command_7

#####		organize file		#####
#echo "mv ./$sample_type/file_to_remove"

