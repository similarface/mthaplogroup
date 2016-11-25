#coding:utf-8
__author__ = 'similarface'

from collections import defaultdict
import math
import json
def getSnpWeightDict():
    '''
    获取权重字典
    算法：F=max(v) F标示某个snp最大出现的频次
    for i in I put:w(i)=10-9lnf(i)/ln(F)
    :param re:
    :param maxcount:
    :return:
    '''
    re=defaultdict(int)
    mttree=json.loads(open('./mttreedb.json','r').read())
    scanmttreesetsnpcount(mttree,re)
    weightdict={}
    maxcount=getMaxItem(re)
    for k,v in re.items():
        #weightdict[k]=10-round(9*math.log10(v)/math.log10(maxcount),1)
        weightdict[k]=10-round(9*math.log(v,math.e)/math.log(maxcount,math.e),1)
    return weightdict


def scanmttreesetsnpcount(tree,re):
    '''
    遍历字典获取snp出现的频率数
    :param tree:树
    :param re:snp 频次字典
    :return:
    '''
    rsid=tree['rsid']
    for rs in rsid:
        re[rs]=re[rs]+1
    children=tree['children']
    if len(children)!=0:
        for child in children:
            scanmttreesetsnpcount(child,re)

def getMaxItem(re):
    '''
    获取最大的次数
    :param re:
    :return:
    '''
    maxcount=0
    for k,v in re.items():
        if v>maxcount:
            maxcount=v
    return maxcount

def getTotalWeightTree(sdata,weightdict,mttree):
    '''
    获取某一个sample所有突变点的权重之和

    :param sdata: 芯片数据
    :param weightdict: 权重字典
    :param mttree:
    :return:所有突变点的权重之和
    '''
    weightTree=0
    rsidset=[]
    getTotalRsidInMtTree(mttree,rsidset)
    for rid in list(set(rsidset)):
        try:
            if rid[-1] == sdata[rid[1:-1]].genotype[0]:
                weightTree=weightTree+weightdict[rid]
        except Exception,e:
            pass
    return weightTree


def getTotalRsidInMtTree(mttree,rsidlist):
    '''
    遍历字典树 获取所有的rsid
    :param mttree:
    :param rsidlist:
    :return:
    '''
    rsid=mttree['rsid']
    if len(rsid)!=0:
        for rs in rsid:
            rsidlist.append(rs)
    children=mttree['children']
    if len(children)!=0:
        for child in children:
            getTotalRsidInMtTree(child,rsidlist)

if __name__=='__main__':
    print(getSnpWeightDict())