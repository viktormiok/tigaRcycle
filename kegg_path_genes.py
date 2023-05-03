from bioservices.kegg import KEGG
import re
import pandas as pd
dirr = '~/Documents/consultation/ChunXia/Rat_timedata/circadian_analysis/input/'

path = pd.read_excel(dirr + 'target_genes.xlsx',
                     sheet_name='KEGG'
       )['KEGG']

def substring(whole, sub1, sub2):
    return whole[whole.index(sub1) : whole.index(sub2)]

def kegg_genes(pathways):
    k = KEGG()
    k.organism = "rno"
    result = {}
    for i in path:
        data = k.get(i[0:8])
        dict_data = k.parse(data)
        if(type(dict_data) == str):
            if("COMPOUND" in dict_data):
                dict_data = substring(dict_data, 
                                      "GENE", 
                                      "COMPOUND"
                )
            else: 
                dict_data = substring(dict_data,
                                      "GENE",
                                      "REFERENCE"
                )
            subStr = re.findall(r' (.+?);',
                                dict_data
            )
            subStr = [x.split() for x in subStr]
        if(type(dict_data) == dict):
            subStr = [x.split() for x in list(dict_data['GENE'])]
        pathgen = []
        for x in range(0, len(subStr)):
            pathgen.append(subStr[x][1])
        result[i] = pathgen
    return result

mydict = kegg_genes(pathways=path)


max_len = []
for key in mydict:
    max_len.append(len(mydict[key]))

for key in mydict:
    la = len(mydict[key])
    if not max_len == la:
        mydict[key].extend(['']*(max(max_len)-la))

df = pd.DataFrame.from_dict(mydict)
df.to_excel(dirr + 'kegg_genes.xlsx',
            header=True,
            index=False
)
