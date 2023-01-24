import os,sys

cohortfile = open(sys.argv[1],'r')
DATE = sys.argv[2]


lines = cohortfile.readlines()

header = open('pre-made-header.vcf','r')
headerlines = header.readlines()

settings = {
    'version': '3.0',
    'stageinmode': 'symlink',
    'population': 'StOlavTAR',
    'annotation': {
        'allele_count': 'AC_StOlavTAR',
        'count_hom': 'Hom_StOlavTAR',
        'count_het': 'Het_StOlavTAR',
        'allele_number': 'AN_StOlavTAR',
        'allele_frequency': 'AF_StOlavTAR',
        'filter': 'filter_StOlavTAR',
        'ALL': 'indications_StOlavTAR'
    },
    'OUTPUT': DATE+'/inhousedb.vcf',
    'OUTPUT_SENSITIVE': DATE+'/inhousedb_sensitive.vcf',
    'INDICATIONS_THRESHOLD': 5
}


def _parse_line(line):
    parts = line.strip().split('\t')
    genotype_field = parts[9:]
    #indications = defaultdict(int)
    count_hom = 0
    count_het = 0
    NoGTlist = ['2','3','4','5','6']
    rmSamples = 0
    
    for genotype in genotype_field:
        gt0, gt1 = genotype.split(':')[0].split('/')
        ## remove multiallec sites.
        if gt0 in NoGTlist or gt1 in NoGTlist:
            rmSamples +=1
            print(gt0,gt1,genotype )
            return None
            #continue
        else:
            tr = {'.': 0, '1': 1, '0': 0}
            gt_sum = tr[gt0] + tr[gt1]
            assert gt_sum in [0, 1, 2]
            #if gt_sum:
            #    indications[gene_panel] += 1
            if gt_sum == 2:
                count_hom += 1
            elif gt_sum == 1:
                count_het += 1
            #else:
            #    print("what")
            #    continue
    actualSamples = len(genotype_field)-rmSamples
    #print(actualSamples)
    allele_count = 2*count_hom + count_het
    allele_number = 2*actualSamples
    allele_frequency = allele_count / float(allele_number)
    if allele_count == 0:
        return None

    db = {
        'CHROM': parts[0].replace("chr",""),
        'POS': parts[1],
        'ID': parts[2],
        'REF': parts[3],
        'ALT': parts[4],
        'QUAL': parts[5],
        'FILTER': "PASS",
        'allele_count': allele_count,
        'count_hom': count_hom,
        'count_het': count_het,
        'allele_number': allele_number,
        'allele_frequency': allele_frequency
    }

    return db

def _create_line(settings,db):
    if db==None:
        return ""
    annotationSet = settings['annotation']
    annotation = {}
    
    annotation['allele_count'] = str(db['allele_count'])
    annotation['count_hom']= str(db['count_hom'])
    annotation['count_het']= str(db['count_het'])
    annotation['allele_number']= str(db['allele_number'])
    annotation['allele_frequency']= '{:.4f}'.format(db['allele_frequency'])
    annotation['filter']= ','.join(db['FILTER'].split(';'))
    annotation['ALL']="ALL:"+str(int(annotation['allele_number'])/2.0)
    db['INFO'] = ""
    #print(annotationSet.keys())
    for key in annotationSet:

        db['INFO'] +=annotationSet[key]+"="+annotation[key]+";"
    db['INFO'] = db['INFO'][0:-1]
    #db['INFO'] = ';'.join(['{}={}'.format(k, v) for k, v in info.iteritems()])
    #db['INFO'] = ';'.join(['{}={}'.format(k, v) for k, v in info])
    db['FORMAT'] = 'GT'
    db['DUMMY'] = './.'

    line = [db[key] for key in ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'DUMMY']]
    line =  '\t'.join(line) + '\n'
    return line



outfile = open(settings['OUTPUT'],'w')
outfile.write("".join(str(s) for s in headerlines))

for line in lines:
    if line.find("#")==0:
        continue
    else:
        db = _parse_line(line)
        line = _create_line(settings,db)
        if line=="":
            continue
        outfile.write(line.strip()+os.linesep)

