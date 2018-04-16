import sys
import urllib.request
import re

files = {
'Athaliana':['ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR292/SRR2922654/SRR2922654.sra'],
'Athaliana-2':['ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR342/SRR342382/SRR342382.sra'],
'Athaliana-epiRIL-01':['ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR292/SRR2922656/SRR2922656.sra'],
'Athaliana-epiRIL-12':['ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR292/SRR2922657/SRR2922657.sra'],
'Athaliana-epiRIL-28':['ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR292/SRR2922658/SRR2922658.sra'],
'Atrichopoda':['ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR116/SRR1168334/SRR1168334.sra'],
'Bdistachyon':['ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR328/SRR3286247/SRR3286247.sra',
'ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR328/SRR3286248/SRR3286248.sra',
'ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR328/SRR3286249/SRR3286249.sra',
'ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR328/SRR3286250/SRR3286250.sra',
'ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR328/SRR3286251/SRR3286251.sra',
'ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR328/SRR3286252/SRR3286252.sra'],
'Eguineensis':['ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR199/SRR1999388/SRR1999388.sra',
'ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR199/SRR1999389/SRR1999389.sra'],
'Esalsugineum':['ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR292/SRR2922631/SRR2922631.sra',
'ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR292/SRR2922632/SRR2922632.sra',
'ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR292/SRR2922633/SRR2922633.sra',
'ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR292/SRR2922634/SRR2922634.sra',
'ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR292/SRR2922635/SRR2922635.sra',
'ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR292/SRR2922636/SRR2922636.sra'],
'Macuminata':['ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR632/SRR6328804/SRR6328804.sra'],
'Osativa':[''],
'Sbicolor':['ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR328/SRR3286309/SRR3286309.sra'],
'Sitalica':['ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR574/SRR5748762/SRR5748762.sra'],
'Sitalica-2':['ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR574/SRR5748763/SRR5748763.sra'],
'Zmays':['ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR207/SRR2079447/SRR2079447.sra']
}

for i in files.get(sys.argv[1]):
	print(re.sub('ftp.*\/','',i))
	urllib.request.urlretrieve(i,re.sub('ftp.*\/','',i))

