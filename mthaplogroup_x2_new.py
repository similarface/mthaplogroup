#!/usr/bin/env python
# coding:utf-8
__author__ = 'similarface'
import snp
import json
from collections import deque
import logging
import logging.handlers
import sys


from WeightCalculation import getSnpWeightDict,getTotalWeightTree

'''
读取23andme的芯片数据
'''

#日子模块
LOG_FILE = 'mtaplogroup.log'
handler = logging.handlers.RotatingFileHandler(LOG_FILE, maxBytes=10240 * 10240, backupCount=1)  # 实例化handler
fmt = '%(asctime)s - %(filename)s:%(lineno)s - %(name)s - %(message)s'
formatter = logging.Formatter(fmt)  # 实例化formatter
handler.setFormatter(formatter)  # 为handler添加formatter
loggername = "mthaplogfile_x2"
logger = logging.getLogger(loggername)  # 获取名为tst的logger
logger.addHandler(handler)  # 为logger添加handler
logger.setLevel(logging.DEBUG)


#mtDNA phylogenetic tree
#加载mtDNA phylogenetic tree json为字典
with open('./mttreedb.json', 'r') as mttreedata:
    mttree = json.loads(mttreedata.read())

#snp的频次字典
weightdict=getSnpWeightDict()


def getMthaologroup(mttree, sdata, result,weightTree):
    '''
    :param ytree: Y染色体的haplogroup树
    :param ysnp: Y染色体突变字典
    :param sdata: 芯片数据
    :param result: 结果
    :return:
    '''
    id = mttree['id']
    pid = mttree['pid']
    name = mttree['name']
    children = mttree['children']
    rsids = mttree['rsid']
    # 匹配的rsid
    matchrsid = []
    nomatchflag=True
    # 不匹配的rsid
    nomathrsid = []
    qfeziweight=0
    qfenmuweight=0
    unkown=0
    unkownflag=True
    unkownweight=0
    for rsid in rsids:
        try:
            if rsid[-1] == sdata[rsid[1:-1]].genotype[0]:
                unkownflag=False
                logger.debug("pos {0},mutation {1},name {2} ,match {3}".format(rsid,sdata[rsid[1:-1]].genotype[0], name,'YES'))
                matchrsid.append(rsid)
                qfeziweight=weightdict[rsid]+qfeziweight
                qfenmuweight=weightdict[rsid]+qfenmuweight
                nomatchflag=False
            else:
                nomatchflag=False
                unkownflag=False
                logger.debug("pos {0},mutation {1},name {2} ,match {3}".format(rsid,sdata[rsid[1:-1]].genotype[0], name,'NO'))
                #weight=weight-1
                qfenmuweight=weightdict[rsid]+qfenmuweight
                nomathrsid.append(rsid)
        except KeyError, e:
            if unkownflag:
                unkown=unkown+1
                unkownweight=(weightdict[rsid])+unkownweight
            qfenmuweight=weightdict[rsid]+qfenmuweight
            nomathrsid.append(rsid)
            try:
                logger.debug("pos {0},mutation {1},name {2} ,match {3}".format(rsid, '?', name,
                                                                                       'UNKOWNU'))
            except KeyError, e:
                try:
                    logger.debug(
                        "pos {0},mutation {1},name {2} ,match {3}".format(rsid, sdata[rsid[1:-1]].genotype[0],
                                                                                  name, 'UNKOWNU'))
                except KeyError, e:
                    logger.debug(
                        "pos {0},mutation {1},name {2} ,match {3}".format(rsid, '?', name, 'UNKOWNU'))
        if not unkownflag:
            unkownweight=0
    if nomatchflag:
        currentfenmuweight=result[pid]['qfenmuweight']
    else:
        currentfenmuweight=result[pid]['qfenmuweight']+qfenmuweight
    #currentfenziweight=result[pid]['qfeziweight']+qfeziweight+unkownweight/pow(2,unkown)
    currentfenziweight=result[pid]['qfeziweight']+qfeziweight
    #currentfenmuweight=result[pid]['qfenmuweight']+qfenmuweight
    #print(id,currentfenziweight,currentfenmuweight)
    if currentfenmuweight!=0:
        result[id] = {'pid': pid, 'name': name, 'matchrsid': matchrsid, 'nomathrsid': nomathrsid,
                  'weight': (0.5*currentfenziweight/currentfenmuweight+2*currentfenziweight/weightTree),
                       "qfeziweight":currentfenziweight,
                  "qfenmuweight":currentfenmuweight}
    else:
        result[id] = {'pid': pid, 'name': name, 'matchrsid': matchrsid, 'nomathrsid': nomathrsid,
                  'weight': 0+2*currentfenziweight/weightTree,
                       "qfeziweight":currentfenziweight,
                  "qfenmuweight":currentfenmuweight}
    #print(id,result[id]['weight'],currentfenziweight,currentfenmuweight,name,weightTree)
    if len(children) != 0:
        for itemdict in children:
            getMthaologroup(itemdict, sdata, result,weightTree)
    else:
        pass

def getPath(mttree, deques, id):
    '''
    根据ID 获取树的路径
    :param y-tree等级树
    :param deques: 是个双端队列 可以前面插入值
    :param id: 树中的某个节点
    :return:返回队列 deques
    '''
    # id=0表示到达树的顶层
    if id == "0":
        return ''
    # 左侧加入该节点的父亲节点
    try:
        deques.appendleft({mttree[id]['name']:{'matchrsid':mttree[id]['matchrsid'],'nomatchrsid':mttree[id]['nomathrsid']}})
    except Exception,e:
        print(id)
    # 获取ytree树的父节点的PID
    pid = mttree[id]['pid']
    getPath(mttree, deques, pid)


def getMaxWeightNode(re):
    # 权重值
    imax = -1000
    # 字典的key
    maxid = 0
    lre={}
    for k,v in re.items():
        lre[k]=v['weight']
    return sorted(lre.iteritems(),key=lambda asd:asd[1],reverse=True)


def getMthaologroupResult(sdata,weightTree):
    # 存储运算结果
    result = {"0": {'matchrsid': None, 'nomathrsid': None, 'weight': 0, 'name': None,'lain':0,'qfeziweight':0,'qfenmuweight':0}}
    getMthaologroup(mttree, sdata, result,weightTree)
    return result


def getMthaologroupPath(sdata,weightTree):
    result = getMthaologroupResult(sdata,weightTree)
    # 获取权重最大的节点
    alltreepath=[]
    i=0
    oldmaxid=-100000
    lenzeropath=-100000
    retreepath=None
    for node in getMaxWeightNode(result):
        i=i+1
        if i>500 :
            break
        maxid, imax = node
        #print(oldmaxid,imax)
        # if oldmaxid>imax:
        #     break
        logger.info('maxid {0} ,imax {1}'.format(maxid, imax))
        treepath = deque([])
        # 获取权重最大的节点的路径
        #getPath(result, treepath, maxid)
        if len(result[maxid]['matchrsid'])!=0:
            print("权重:"+str(imax))
            logger.info('maxid {0} ,imax {1}'.format(maxid, imax))
            oldmaxid=imax
            getPath(result, treepath, maxid)
            #打印匹配信息
            for item in treepath:
                for k,v in item.items():
                    print(k,v['matchrsid'],v['nomatchrsid'])

            logger.info('->'.join(getFormatYhaologroupPath(treepath)))
            #logger.info('>'.join(getFormatYhaologroupPath(treepath)[1]))
            zeropath=getFormatYhaologroupPath(treepath)[1]
            #print '->'.join(getFormatYhaologroupPath(treepath)[0])
            #print '->'.join(getFormatYhaologroupPath(treepath)[0])
            if lenzeropath<len(zeropath):
                lenzeropath=len(zeropath)
                retreepath=treepath

            #print(getZeroCount(zeropath),getSumLain(zeropath),len(zeropath))
            #print(treepath)
            alltreepath.append(treepath)

            print '->'.join(getFormatYhaologroupPath(treepath))
            print("x2\t"+getFormatYhaologroupPath(treepath)[-1])
    return retreepath


def getFormatYhaologroupPath(treepath):
    re = []
    for item in treepath:
        for k,v in item.items():
            re.append(k)
    return re

def main_x2(filename):
    #更新位置信息
    #update23andmeFilePos(filename)
    sdata = snp.reader_dict_mtpos(filename)
    #所有命中snp的权重和
    weightTree=getTotalWeightTree(sdata,weightdict,mttree)
    if weightTree>0:
        getMthaologroupPath(sdata,weightTree)
    else:
        print('请检查文件')

if __name__ == "__main__":
    #filename='/Users/similarface/Documents/user2533_file1566_yearofbirth_1969_sex_XX.23andme.txt'
    filename='/Users/similarface/Documents/FangGE.txt'

    if len(sys.argv)!=1:
        print("param error!","Use: python mthaplogroup andmefile")
    else:
        main_x2(filename)