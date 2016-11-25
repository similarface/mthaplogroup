
import csv
import MySQLdb
reader = csv.reader(file('/Users/similarface/Documents/batch2.csv/batch22.csv', 'rb'))

cursor=MySQLdb.connect(host='192.168.30.252',port=3306,db='gendb',user='dna',passwd='dna',charset='utf8').cursor()

for line in reader:
    print line
    SQL="insert into t_20160107_batch_2_stable(flagt   ,chrpos   ,chr      ,bpos     ,epos     ,ref      ,alt      ,rsid     ,dbsource ,ntype     ,cname    ,rsids  ) values(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)"
    cursor.execute(SQL, line)

cursor.close()
