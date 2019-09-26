import sys,os
import pandas as pd
A = pd.read_table(sys.argv[1], header=None, sep='\s+',chunksize =50,lineterminator='\n')
A_df = pd.concat(A, ignore_index=True)
A_df.columns = ['chr', 'start','stop','reads','temp']
A_df['id'] = A_df[['chr', 'start','stop']].apply(lambda x: '_'.join(str(value) for value in x), axis=1)
B = pd.read_table(sys.argv[2], header=None, sep='\s+',chunksize =50,lineterminator='\n')
B_df = pd.concat(B, ignore_index=True)
B_df.columns = ['chr', 'start','stop','reads','temp']
B_df['id'] = B_df[['chr', 'start','stop']].apply(lambda x: '_'.join(str(value) for value in x), axis=1)
PAV=sys.argv[3]
output="/cluster/project/gdc/people/crimjhim/"+PAV
#output="/Users/rimjhim/Documents/WORK/Arabis_resequencing/test_zygosity/new/"+PAV
#print output
#print output
os.chdir(output)
for filename in os.listdir(output):
	if filename.endswith(".txt"):
		#print filename
		locus=filename.split(".txt")[0]
		C_df = pd.read_table(filename, header=None, sep='\s+')
		C_df.columns = ['chr', 'start','stop']
		C_df['id'] = C_df[['chr', 'start','stop']].apply(lambda x: '_'.join(str(value) for value in x), axis=1)
		new_df = pd.merge(A_df, C_df, left_on=['id'], right_on = ['id'])
		new2_df = pd.merge(B_df, C_df, left_on=['id'], right_on = ['id'])
		
		if not new_df.empty:
			readlist=str(','.join(map(str,new_df['reads'])))
			myList = [i[:-2] for i in readlist.split (',')] 
		
			slen=str(len(set(myList)))
		else:
			slen=0
		if not new2_df.empty:
			readlist2=str(','.join(map(str,new2_df['reads'])))
			myList2 = [i[:-2] for i in readlist2.split (',')] 
		
			nslen=str(len(set(myList2)))
		else:
			nslen=0
		if (int(slen)+int(nslen))!=0:
			zyg=float(slen)/(float(slen)+float(nslen))
		else:
			zyg=0
		print locus+'\t'+str(slen)+'\t'+str(nslen)+'\t'+str(zyg)

