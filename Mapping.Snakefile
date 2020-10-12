import os                                                                       
import glob                                                                     
import re                                                                       
import paramiko                                                                 
                                                                                
import pandas as pd                                                             
import minus80 as m80                                                           
                                                                                
from urllib.parse import urlparse                                               
                                                                                
                                                                                
# ----------------------------------------------------------                    
#       Set up Minus80 Objects                                                  
# ----------------------------------------------------------                    
MFCohort = m80.Cohort('MuscleFat')                                              
SAMPLES = MFCohort.names                                                        
                                                                                
# ----------------------------------------------------------                    
#       One ring to Rule them all                                               
# ----------------------------------------------------------                    
                                                                                
rule all:                                                                       
    input:                                                                      
        'data/SalmonQuant/MFCohort.tsv'                                         
        #expand('data/SalmonQuant/{sample}',sample=SAMPLES)                     
                                                                                
# ----------------------------------------------------------                    
#       Rules                                                                   
# ----------------------------------------------------------                    
                                                                                
rule trim_reads:                                                                
    input:                                                                      
        R1 = lambda wc: sorted([urlparse(x).path for x in MFCohort[wc.sample].files if 'R1' in x]),
        R2 = lambda wc: sorted([urlparse(x).path for x in MFCohort[wc.sample].files if 'R2' in x])
    output:                                                                     
        R1 = 'data/trimmedfastq/{sample}.R1.fastq.gz',                          
        R2 = 'data/trimmedfastq/{sample}.R2.fastq.gz'                           
    run:                                                                        
        shell('''                                                               
            AdapterRemoval \                                                    
            --file1 {input.R1} \                                                
            --file2 {input.R2} \                                                
            --combined-output \                                                 
            --output1 {output.R1} \                                             
            --output2 {output.R2} \                                             
            --gzip \                                                            
            --trimns \                                                          
            --trimqualities \                                                   
            --minquality 10 \                                                   
        ''')                                                                    
                                                                                
rule SALMON_mapping:                                                            
    input:                                                                      
        R1 = 'data/trimmedfastq/{sample}.R1.fastq.gz',                          
        R2 = 'data/trimmedfastq/{sample}.R2.fastq.gz'                           
    output:                                                                     
        directory('data/SalmonQuant/{sample}'),                                 
        'data/SalmonQuant/{sample}/quant.sf'                                    
    run:                                                                        
        shell('''                                                               
            salmon quant \                                                      
            -i data/RefGen/SALMON_INDEX \                                       
            -l A \                                                              
            -1 {input.R1} \                                                     
            -2 {input.R2} \                                                     
            --validateMappings \                                                
            -o {output}                                                         
        ''')                                                                    
                                                                                
rule make_TPM_matrix:                                                           
    input:                                                                      
        quant_files = expand('data/SalmonQuant/{sample}/quant.sf',sample=SAMPLES)
    output:                                                                     
        tpm_matrix = 'data/SalmonQuant/MFCohort_TPM.tsv',                       
        fat_matrix = 'data/SalmonQuant/Fat_TPM.tsv',                            
        muscle_matrix = 'data/SalmonQuant/Muscle_TPM.tsv'                       
    run:                                                                        
        tpm = None                                                              
        for n in input.quant_files:                                             
            name = n.split('/')[2]                                              
            #x = pd.read_table(f"data/SalmonQuant/{n}/quant.sf")                                         
            x = pd.read_table(n)                                                
            x = x.loc[:,['Name','TPM']]                                         
            x.rename(columns={'Name':'gene','TPM':name},inplace=True)           
            x.set_index('gene',inplace=True)                                    
            if tpm is None:                                                     
                tpm = x                                                         
            else:                                                               
                tpm = tpm.join(x)                                               
        # Create the whole matrix                                               
        tpm.to_csv(output.tpm_matrix,sep='\t')                                  
        # Create the Fat & Muscle Matrix                                        
        tpm.loc[:,[x.endswith('F') for x in tpm.columns]].to_csv(output.fat_matrix,sep='\t')
        tpm.loc[:,[x.endswith('M') for x in tpm.columns]].to_csv(output.muscle_matrix,sep='\t')
