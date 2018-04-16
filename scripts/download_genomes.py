import sys
from subprocess import call

#Determine source database
if sys.argv[2] == "phytozome":
    #Create cookies
    url1='https://signon.jgi.doe.gov/signon/create'
    login='login=<USERNAME>'
    password='password=<PASSWORD>'
    call(['curl',url1,'--data-urlencode',login,'--data-urlencode',password,'-c','cookies'])
    url2='https://genome.jgi.doe.gov/ext-api/downloads/get_tape_file?blocking=true&url=/PhytozomeV12/download/_JAMO/'
elif sys.argv[2] == "ensemble":
    url2="ftp://ftp.ensemblgenomes.org/pub/plants/release-37/"
else:
    print("A source site is needed!")

#List files
genome={
'Athaliana.fa':'587b0adf7ded5e4229d885ab/Athaliana_447_TAIR10.fa.gz',
'Atrichopoda.fa':'56981cdc0d87851ee9727d1f/Atrichopoda_291_v1.0.fa.gz',
'Bdistachyon.fa':'597f9cc67ded5e0452b3f390/SbicolorRio_468_v2.0.fa.gz',
'Eguineensis.fa':'',
'Esalsugineum.fa':'56981cdc0d87851ee9727d15/Esalsugineum_173_v1.fa.gz',
'Osativa.fa':'',
'Sbicolor.fa':'',
'Zmays':'fasta/zea_mays/dna/Zea_mays.AGPv4.dna.toplevel.fa.gz'
}

gff={
'Athaliana.gff':'5854743b7ded5e78cff8c475/Athaliana_167_TAIR10.gene_exons.gff3.gz',
'Atrichopoda.gff':'56981cdc0d87851ee9727d18/Atrichopoda_291_v1.0.gene_exons.gff3.gz',
'Bdistachyon.gff':'597f9cc57ded5e0452b3f38e/SbicolorRio_468_v2.1.gene_exons.gff3.gz',
'Eguineensis.gff':'',
'Esalsugineum.gff':'56981cdc0d87851ee9727d0e/Esalsugineum_173_v1.0.gene_exons.gff3.gz',
'Osativa.gff':'',
'Sbicolor.gff':'',
'Zmays.gff':'gff3/zea_mays/Zea_mays.AGPv4.37.gff3.gz'
}

#Download files
if genome.get(sys.argv[1]+'.fa'):
    f=open(sys.argv[1]+'.fa.gz','w')
    call(['curl',url2+genome.get(sys.argv[1]+'.fa'),'-b','cookies'],stdout=f)
    f.close()
if gff.get(sys.argv[1]+'.gff'):
    f=open(sys.argv[1]+'.gff.gz','w')
    call(['curl',url2+gff.get(sys.argv[1]+'.gff'),'-b','cookies'],stdout=f)
    f.close()
