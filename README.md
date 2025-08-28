Code used for analyses in Amourda et al, Myc sustains Regional and Sex-biased Organ Zonation in the Drosophila Intestine

Lipid droplet counting pipeline:
Pipeline for identifying lipid droplets and measuring number, volume and intensity of each droplet and assigning this along gut length. Also measures overall average intensity of fluorescence signal along gut length.

System requirements:
Scripts written and tested on an Apple MacBook Pro 2023 with an Apple M3 Pro Chip and 18GB of RAM and a Windows virtual machine with 2.8GHz Intel(R) Xeon(R) Gold CPU and 64GB RAM using FIJI version 2.14.0/1.54p and R version 4.4.0.

Installation: 
This analytical code does not require installation beyond standard installation of R packages as specified in the script and requires FIJI plugins: 3D image suite v4.1.7 and MorphoLibJ v1.6.4.

Input:
Confocal .lif file containing stacks of gut rolls imaged with high resolution Leica 'Lightning' mode

To run:
First step - run 'LDpipeline_Run_all.ijm' macro in FIJI and follow instructions provided by macro
Second step - run LDpipline R scripts parts 1-5 in R, changing directories as specified in scripts

Expected output:
Graphs of lipid droplet number, average droplet volume, average droplet intensity along gut length, graphs of average intensity along gut length and output of statistical tests on area under the curve for the R2 and R4 gut regions.
