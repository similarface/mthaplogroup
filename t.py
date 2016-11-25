__author__ = 'similarface'
import json
with open('./mttreedb.json', 'r') as mttreedata:
    mttree = json.loads(mttreedata.read())

def jsonprint(it):
    if len(it['rsid'])!=0:
        rsids=it['rsid']
        #print rsids
        for r in rsids:
            res=list(r)
            pre=""
            pos=""
            end=""
            flag=True
            for c in res:
                if c.isdigit():
                    pos=pos+c
                    flag=False
                else:
                    if flag:
                        pre=pre+c
                    else:
                        end=end+c
            print '\t'.join(['MT',pos,pre,'->',end])


    if it['children']!=None:
        for x in it['children']:
            jsonprint(x)

if __name__ == "__main__":
    jsonprint(mttree)