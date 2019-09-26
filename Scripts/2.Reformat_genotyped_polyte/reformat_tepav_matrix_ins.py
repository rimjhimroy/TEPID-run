#! /usr/bin/env python

def convert_matrix(infile, names,refdict):
    """
    convert TEPAV genotyped file to presence / absence matrix
    """
    
    if infile.endswith(".gz"):
        inf = gzip.open(infile, "rb")
    
    else:
        inf = open(infile, "r")
    for line in inf:
        line = line.rsplit()
        if line[0] == "chromosome":
            pass  # header
        else:
            

            chrom = line[0]
            start = line[1]
            stop = line[2]
            TE = line[6]
            TEID=line[3]+'_'+line[4]+'_'+line[5]
            #print TEID
            n_te = str(len(TE.split(",")))
            ntes=TE.split(",")
            #print ntes
            tes=filter(lambda x:TEID in x, ntes)
            if not tes: tes=ntes
            #print tes
            
            tefam=[]
            tetribe=[]
            
            
            for i in xrange(len(tes)):
                
                tefam.append(refdict[tes[i]][0])
                
                tetribe.append(refdict[tes[i]][1])
            
            
            tribe=list(set(tetribe))
            if 'Unknown' in tribe:
                tribe.remove('Unknown')
            if not tribe:
                tribe.append('Unknown')
            
            pos = line[7].split(",")
            neg = line[8].split(",")
            te_id = "\t".join([chrom, start, stop])
            status = get_status(pos, neg, names)
            column_ordered = []
            for i in names:
                column_ordered.append(status[i])
            noNA = filter(lambda x: x != "NA", status.values())            
            noNA = map(int, noNA)
            pos_count = sum(noNA)
            l = len(noNA)
            neg_count = l - pos_count
            TE_present=pos_count
            TE_absent=neg_count
            if(pos_count < neg_count):
                Minor_allele="presence"

            else:
                Minor_allele="absence"
#print Minor_allele
            if (TE_absent < 61):
                Absence_classification="True deletion"
            elif (TE_absent > 244):
                Absence_classification="No insertion"
            else:
                Absence_classification="NA"
            original_call_deletion = 'F'
            MAF=float(min(TE_present, TE_absent))/304
            #print int(min(TE_present, TE_absent)) ,MAF
            if(MAF < 0.025):
                Frequency_classification = "Rare"
            else:Frequency_classification ="Common"
            print(te_id + "\t" + TE + "\t" + ",".join(tefam) + "\t" +",".join(tribe) + "\t" +n_te + "\t"  + str(pos_count) + "\t" + str(neg_count) + "\t" + str(Minor_allele) + "\t" +original_call_deletion + "\t" +str(Absence_classification) + "\t" +str(MAF) + "\t" +str(Frequency_classification) + "\t"+"\t".join(column_ordered))
    inf.close()


def create_names_dict(infile):
    """
    read accession names into list
    """
    return [name.strip("\n") for name in open(infile, "r")]


def get_status(pos, neg, names):
    """
    create dictionary of TEPAV status for each accession (1/0/NA)
    """
    status = {}
    for i in names:
        #print str(i) +'\n'+ str(pos) +'\n'+ str(neg)+'\n'+'\n'
        if i in pos:
            status[i] = "1"
        elif i in neg:
            status[i] = "0"
        else:
            status[i] = "NA"
    return status


if __name__ == "__main__":
    from argparse import ArgumentParser
    import gzip
    import re,csv

    parser = ArgumentParser(description='Convert TEPAV genotyped file to presence/absence matrix')
    parser.add_argument('-n', '--names', help='All sample names newline separated', required=True)
    parser.add_argument('-i', '--infile', help="Input TEPAV file", required=True)
    options = parser.parse_args()

    names = create_names_dict(options.names)
    print "chrom" +"\t" + "start" + "\t" + "stop" + "\t" + "TEID" + "\t" + "FAMILY" + "\t" +"TRIBE" + "\t" +"number_of_te" + "\t" + "TE_present" + "\t" + "TE_absent" + "\t" + "Minor_alleles" + "\t" + "original_call_deletion" + "\t" +"Absence_classification" + "\t" +"MAF" +"\t" +"Frequency_classification"+"\t" +"\t".join(names)
    ref = list(csv.reader(open("reference.full.solo.repeat", 'rb'), delimiter='\t'))
    refdict=dict()
    for i in xrange(len(ref)):
        
        key = str(ref[i][1].strip())
        #print key
        family = ref[i][0].strip()
        tribe=ref[i][2].strip()
        if key in refdict :
            refdict[key].append([family,tribe])
        else:
            refdict[key]=([family,tribe])

convert_matrix(options.infile, names,refdict)
