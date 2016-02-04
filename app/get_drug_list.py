from collections import Counter
import numpy as np

def get_drug_list(df):
    pri_list=Counter(df.Primary_site).keys()
    sec_list=['Primary tumor','Metastatic tumor', 'Undetermined']

    mydict=dict()
    for a in pri_list:
        mydict[a]=dict()

    for a in pri_list:
        for b in sec_list:
            mydict[a][b]=np.array(df.loc[(df.Primary_site==a)&(df.Tumour_origin==b),['cell_line_name','Comments']]).tolist()

    return mydict
