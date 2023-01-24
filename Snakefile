#######################
#    Variant frequent database pipeline - 
#
#    Copyright (C) 2023 Ashish Kumar Singh, Jostein Johansen - St.Olavs hospital - Bioinformatics core facility - NTNU - Norway
#
#    The Variant frequent database pipeline  and all scripts available is free: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This pipeline is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
########################

configfile: "config.yaml"

import pandas as pd
import io
from datetime import date


HOME = config['HOME']
REFERENCE = config['REFERENCE']
SCRIPTS = config['SCRIPTS']

EXPORT = config['EXPORT']

GATK = config['GATK']
RESULTS = config['RESULTS']
DATABASE = config['DATABASE']



PREVDATE = os.listdir("RUNS")[-1]



today = date.today()
DATE = today.strftime("%Y-%m-%d")
if DATE!=PREVDATE:
    WORKDIR = HOME+"/RUNS/"+str(DATE)
    print("Starts new analysis" + DATE)
    os.system("mkdir -p "+WORKDIR)
else:
    DATE=PREVDATE
    WORKDIR = HOME+"/RUNS/"+str(PREVDATE)


def checkSamples(RESULTS,DATABASE,WORKDIR):
    gvcfsamples = os.listdir(DATABASE)
    checkfiles = os.listdir(RESULTS)
    newcheckfiles = []
    for folder in checkfiles:
        if folder.find("TAR")==-1:
            continue
        if folder.find("NGS-Quality-TAR")!=-1:
            continue
        try:
            samples = os.listdir(RESULTS+"/"+folder+"/results/samples/")
            for s in samples:
                vcffolderfiles = os.listdir(RESULTS+"/"+folder+"/results/samples/"+s+"/vcf")
                for v in vcffolderfiles:
                    if v.find("_raw_SNP_INDEL.g.vcf")!=-1:
                        if v in gvcfsamples:
                            continue
                        else:
                            if s.find("100016253159")==-1:
                                newcheckfiles.append(RESULTS+"/"+folder+"/results/samples/"+s+"/vcf/"+v)
        except:
            "OLD FOLDER STRUCTURE"
    logfile = open('LOG/'+DATE+'-log.txt','w')
    if len(newcheckfiles)>0:
        for n in newcheckfiles:
            logfile.write("added sample: "+n+"\n")
            os.system("cp "+n+" "+DATABASE)
        samples = os.listdir(DATABASE)
        allGVCFlist = open(WORKDIR+"/ALL_gVCF.list",'w')
        for s in samples:
            if s.find(".idx")==-1:
                allGVCFlist.write(DATABASE+"/"+s+"\n")
        allGVCFlist.close()
    return samples


samples = checkSamples(RESULTS,DATABASE,WORKDIR)


rule final_results:
    input:
        #WORKDIR+"/combine_vcfs_done"
        WORKDIR+"/internal_db_done"
        


rule combine_vcfs:
    input:
        WORKDIR
    output:
        touch(WORKDIR+"/combine_vcfs_done")
    threads:
        48
    shell:
        """
        gatk CombineGVCFs \
        --reference {REFERENCE} \
        --variant {input}/ALL_gVCF.list \
        --output  {input}/ALL_gVCF_gatk-CombineGVCFs.g.vcf.gz
        """

rule genotype_gvcfs:
    input:
        WORKDIR+"/combine_vcfs_done"
    output:
        touch(WORKDIR+"/genotype_vcfs_done")
    params:
        WORKDIR
    threads:
        48
    shell:
        """
        gatk GenotypeGVCFs \
        --reference {REFERENCE} \
        --variant {params}/ALL_gVCF_gatk-CombineGVCFs.g.vcf.gz \
        --output {params}/ALL_gVCF_gatk-CombineGVCFs_GenotypeGVCFs.g.vcf.gz \
        --allow-old-rms-mapping-quality-annotation-data 
    """


rule decompose_combined_vcf:
    input:
        WORKDIR+"/genotype_vcfs_done"
    output:
        touch(WORKDIR+"/decompose_combined_vcf_done")
    shell:
        """
       vt decompose -s {WORKDIR}/ALL_gVCF_gatk-CombineGVCFs_GenotypeGVCFs.g.vcf.gz -o {WORKDIR}/ALL_gVCF_gatk-CombineGVCFs_GenotypeGVCFs_vt-decomposed.g.vcf.gz
        """

rule normalize_decomposed_vcf:
    input:
        WORKDIR+"/decompose_combined_vcf_done"
    output:
        touch(WORKDIR+"/normalize_decomposed_vcf_done")
    shell:
        """
        vt normalize -r {REFERENCE} {WORKDIR}/ALL_gVCF_gatk-CombineGVCFs_GenotypeGVCFs_vt-decomposed.g.vcf.gz -o {WORKDIR}/ALL_gVCF_gatk-CombineGVCFs_GenotypeGVCFs_vt-decomposed_normalized.g.vcf.gz
        """

rule index_normalized_vcf:
    input:
        WORKDIR+"/normalize_decomposed_vcf_done"
    output:
        touch(WORKDIR+"/index_normalized_vcf_done")
    shell:
        """
        gatk IndexFeatureFile -I {WORKDIR}/ALL_gVCF_gatk-CombineGVCFs_GenotypeGVCFs_vt-decomposed_normalized.g.vcf.gz
        """

rule indbcreator:
    input:
        WORKDIR+"/index_normalized_vcf_done"
    output:
        touch(WORKDIR+"/internal_db_done")
    shell:
        """
        zcat {WORKDIR}/ALL_gVCF_gatk-CombineGVCFs_GenotypeGVCFs_vt-decomposed_normalized.g.vcf.gz > {WORKDIR}/ALL_gVCF_gatk-CombineGVCFs_GenotypeGVCFs_vt-decomposed_normalized.g.vcf
        python indbcreator_stolav_normalized.py {WORKDIR}/ALL_gVCF_gatk-CombineGVCFs_GenotypeGVCFs_vt-decomposed_normalized.g.vcf {DATE}
        bgzip -c {WORKDIR}/inhousedb.vcf > {WORKDIR}/inhousedb.vcf.gz
        tabix {WORKDIR}/inhousedb.vcf.gz
        mkdir -p EXPORT/inhousedb_{DATE}
        cp {WORKDIR}/inhousedb.vcf.gz* EXPORT/inhousedb_{DATE}/ 
        """

