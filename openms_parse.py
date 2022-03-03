import re,sys
import numpy as np
import csv
from itertools import islice
openms = sys.argv[1]
unimod = sys.argv[2]
                
def string_parse(mod):
    mod1 = re.sub(u"\)*", "", mod).replace("Protein N-term","*n").replace("(", ",").replace(" ","").split(",")
    mod2 = np.array(mod1)
    ptm = []
                
    with open(unimod,newline = '',encoding = 'utf-8') as f:    
        reader = csv.reader(f)
        for row in reader: 
            for i in range(0,len(mod2)):
                if row[0] == mod2[i]:
                    string = row[3] + "," +mod2[i+1]
                    ptm.append(string)
    return ptm


if __name__=="__main__":
    with open(openms,newline = '',encoding = 'utf-8') as f:
            reader = csv.reader(f, delimiter = '\t')
            for row in islice(reader, 1, None):
                fixptm = string_parse(row[2])
                varptm = string_parse(row[3])
                row[2] = ''
                row[3] = ''
                for i in range(0, len(fixptm)):
                    if (i == len(fixptm) - 1):
                        row[2] = row[2] + '--fixed-mod ' + fixptm[i]
                    else:
                        row[2] = row[2] + '--fixed-mod ' + fixptm[i] + ' '
                for j in range(0, len(varptm)):
                    if (j == len(varptm) - 1):
                        row[3] = row[3] + '--var-mod ' + varptm[j]
                    else:
                        row[3] = row[3] + '--var-mod ' + varptm[j] + ' '
                r = open('params.csv', 'a', newline = '')
                writer = csv.writer(r, delimiter = '\t')
                writer.writerow(row)
                r.close()
    f.close()