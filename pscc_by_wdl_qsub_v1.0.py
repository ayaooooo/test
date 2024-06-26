#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
@Author: 满建芬（manjianfen@bgi.com）
@Date:2022-07-21
'''

import sys
import argparse
import os
import logging
import json



rootDir = sys.path[0]

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--work", action="store", dest="work_path", help="out work path")
parser.add_argument("-a", "--auto", action="store_true", dest="auto", help="auto qsub")
args = parser.parse_args()

work_path = args.work_path
auto = args.auto
if work_path is None:
	parser.print_help()
	print('缺少必要参数 -o')
	sys.exit(-1)
else:
	if not os.path.exists(work_path):
		print('流程输出结果目录不存在：' + work_path)
		sys.exit(-1)

	fqpath = os.path.join(work_path, "sequence_info.list")
	infoFile = os.path.join(work_path, "info.txt")
	if not os.path.exists(fqpath) or not os.path.exists(infoFile):
		print("缺少输入文件seque_nce-info.list或者info.txt")
		sys.exit(-1)

logPath = os.path.join(work_path, 'Log')
if not os.path.exists(logPath):
	os.makedirs(logPath)
	logging.basicConfig(filename=os.path.join(logPath, 'log.log'), level=logging.DEBUG, format='%(levelname)s\t%(asctime)s\t%(filename)s\t[line:%(lineno)d]%(message)s', datefmt='%Y-%m-%d %H:%M:%S')


fqpathDict = dict()
infoDict = dict()

#sysId UID slide lane barcode fq1Path fq2Path fq1Stat fq2Stat summaryReport barcodeStat 
##sampleLib seqPlatform barcode fqpath fqpath2
num_sequence_info = 0
with open(fqpath, "r") as f:
	head = f.readline().rstrip().split("\t")
	if len(head) < 5:
		logging.error('sequence_info.list 文件列数不正确：' + fqpath)		
	fqdata = [l.rstrip().split("\t") for l in f.readlines()]
for data in fqdata:
	fqpathDict[data[0]] = [data[4],data[5],data[6]]#sampeID barcode fqpath fqpath2
	num_sequence_info = num_sequence_info + 1
f.close()

#sampleID	sampleLib seqPlatform gender	productClass	hillsID	hospitalName	productName	productID	platform	poolingLib
product = {
'DX0007':"100K",
'DX0006':"1M",
'DX1379':"100K",
'DX1446':"100K",
'DX1445':"1M",
'RD0005':"1M",
'RC0013':"100K",
'SD0283':"100K",
'DX1679':"100K",
'DX1930':"100K",
'RFP006':"100K",
'DX2093':"100K",
'DX1932':"100K",
'DX1932-A':"100K",
'ZK0004':"100K",
'ZK0005':"100K",
'ZK0006':"100K",
'ZK0007':"100K"
}

num_info = 0
with open(infoFile, "r") as ff:
	head = ff.readline().rstrip().split("\t")
	if len(head) < 5:
		logging.error('info.txt 文件列数不正确：' + infoFile)
	infodata = [l.rstrip().split("\t") for l in ff.readlines()]
for data in infodata:
	infoDict[data[0]] = [data[3],product[data[7]],data[2],data[7]]#gender 100K/1M productID seqPlatform
	if data[0] in fqpathDict:
		fqpathDict[data[0]] = [data[1],data[2]] + fqpathDict[data[0]]
	else:
		logging.error('sequence-info.list, info.txt样本名不一致')
	num_info = num_info + 1
ff.close()

#sampleLib seqPlatform barcode fqpath fqpath2

if num_sequence_info != num_info:
	logging.error('sequence-info.list, info.txt 文件行数不一致')

arg_dic = {
"bin_soap":os.path.join(rootDir,"software","SOAP"),
"bin_hg19Index":os.path.join(rootDir,"SE","Hg19_reference","human.fa.index"),
"bin_java":os.path.join(rootDir,"software","miniconda3","envs/kyVirus/bin","java"),
"bin_cromwell":os.path.join(rootDir,"software","cromwell.jar"),
"sc_soap":os.path.join(rootDir,"script","soap.wdl"),
"bin_fqsta":os.path.join(rootDir,"SE","fq.sta.pl"),
"sc_fastq1":os.path.join(rootDir,"script","fastq-1.wdl"),
#"bin_perl":os.path.join(rootDir,"software","miniconda3","bin","perl"),
"bin_perl":"/home/bgiky/env/perl/bin/perl",
"bin_rmdup":os.path.join(rootDir,"SE","rmDup.pl"),
"bin_SOAPcoverage":os.path.join(rootDir,"software","SOAPcoverage"),
"bin_hg19":os.path.join(rootDir,"SE","Hg19_reference","human.fa"),
"bin_Nregion":os.path.join(rootDir,"SE","HG19.N_region"),
"sc_rmdup":os.path.join(rootDir,"script","rmdup.wdl"),
"bin_soapSta":os.path.join(rootDir,"SE","cov.soap.sta.pl"),
"sc_fastq2":os.path.join(rootDir,"script","fastq-2.wdl"),
"bin_DataProduction":os.path.join(rootDir,"SE","Data_production.pl"),
#"bin_pscc":os.path.join(rootDir,"pscc.py"),
"sc_DataProduction":os.path.join(rootDir,"script","Data_Production.wdl"),

"gender_sh":os.path.join(rootDir,"script","gender.sh"),
"product_sh":os.path.join(rootDir,"script","product.sh"),
"bin_hg19":os.path.join(rootDir,"SE","Hg19_reference","human.fa"),
"bin_soap2ext":os.path.join(rootDir,"PSCC","Soap2Ext.pl"),
"bin_tagsgc":os.path.join(rootDir,"PSCC","tags.gc.info.pl"),
"bin_Aneuploid":os.path.join(rootDir,"PSCC","Aneuploid.pl"),
"bin_winRatio":os.path.join(rootDir,"PSCC","winRatio","winRatio_"),
"bin_RatioLoose":os.path.join(rootDir,"PSCC","Ratio.Loose.pl"),
"bin_Segmentation":os.path.join(rootDir,"PSCC","Segmentation.pl"),
"bin_hg19Length":os.path.join(rootDir,"PSCC","Attachment","hg19.length"),
"bin_hg19Nregion":os.path.join(rootDir,"PSCC","Attachment","hg19.N_region"),
"bin_cnvFilter":os.path.join(rootDir,"PSCC","CNVFilter.pl"),
"bin_spValue":"1e-3",
"bin_ppValue":"1e-1",
"bin_Graph":os.path.join(rootDir,"PSCC","Graph.pl"),
"sc_pscc":os.path.join(rootDir,"script","pscc.wdl"),
"bin_zibizslope":os.path.join(rootDir,"SLOPE-code-use-now","bin","zibiz-slope_P_v18.pl"),
"bin_find_merge":os.path.join(rootDir,"SLOPE-code-use-now","bin","find_merge_v2.pl"),
"bin_CNVFilter":os.path.join(rootDir,"SLOPE-code-use-now","bin","CNVFilter.pl"),
"sc_slope":os.path.join(rootDir,"script","slope.wdl"),

"bin_EXT_SPLIT":os.path.join(rootDir,"SLOPE-code-use-now","bin","EXT_SPLIT_individual_v1.pl"),
"bin_slopeWinRatio":os.path.join(rootDir,"SLOPE-code-use-now","control_cnv","winRatio_"),
"bin_GCcorrection":os.path.join(rootDir,"SLOPE-code-use-now","bin","GC_correction.pl"),
"bin_splitRatio":os.path.join(rootDir,"SLOPE-code-use-now","bin","split_for_ratio.pl"),
"sc_GC_correction_ratio":os.path.join(rootDir,"script","GC_correction_ratio.wdl"),
"bin_zibizslope":os.path.join(rootDir,"SLOPE-code-use-now","bin","zibiz-slope_P_v18.pl"),
"bin_find_merge":os.path.join(rootDir,"SLOPE-code-use-now","bin","find_merge_v2.pl"),
"bin_CNVFilter":os.path.join(rootDir,"SLOPE-code-use-now","bin","CNVFilter.pl"),
"sc_slope":os.path.join(rootDir,"script","slope.wdl"),
"bin_python":os.path.join(rootDir,"software","miniconda3","envs/kyCNV/bin","python"),
"bin_gencheck":os.path.join(rootDir,"AIO_sourcecode","cnv","gencheck.py"),
"bin_merge_5M_pro":os.path.join(rootDir,"AIO_sourcecode","result","merge_5M_pro.py"),
"bin_merge100K":os.path.join(rootDir,"AIO_sourcecode","result","merge_100K_v1.py"),
"bin_hg19Nregion":os.path.join(rootDir,"PSCC","Attachment","hg19.N_region"),
"bin_highConfidence100K":os.path.join(rootDir,"AIO_sourcecode","result","high_confidence_100K.py"),
"bin_hg19Length":os.path.join(rootDir,"PSCC","Attachment","hg19.length"),
"bin_hg19cytoBand":os.path.join(rootDir,"PSCC","Attachment","hg19.cytoBand"),
"bin_GraphReport":os.path.join(rootDir,"PSCC/Graph_report.pl"),
"bin_AIOdata":os.path.join(rootDir,"AIO_sourcecode","dat"),
"sc_step4":os.path.join(rootDir,"script","step4.wdl"),
"bin_python2":os.path.join(rootDir,"software","miniconda3","envs/py2/bin","python2"),
"bin_compareResult":os.path.join(rootDir,"AIO_sourcecode","result","compare_result_pro.py"),
"bin_mergepng":os.path.join(rootDir,"AIO_sourcecode","result","pscc_merge_png.py"),
"sc_result":os.path.join(rootDir,"script","result.wdl"),
"bin_Annotation":os.path.join(rootDir,"Annotation_V3.0","KY_Annotation.py"),
"sc_anno":os.path.join(rootDir,"script","Anno.wdl"),
"bin_kneaddata":os.path.join(rootDir,"software/miniconda3/envs/kyVirus/bin/kneaddata"),#/mnt/nas/pipelines/CNVseq100k_wdl/software/miniconda3/envs/kyVirus/bin/bowtie2
"bin_kneaddataDatabase":os.path.join(rootDir,"kyVirus","data","kneaddata_database"),
"bin_trimmomatic":os.path.join(rootDir,"software/miniconda3/envs/kyVirus/share/trimmomatic-0.39-1"),
"bin_software":os.path.join(rootDir,"software/miniconda3/envs/kyVirus/bin"),
"bin_bowtie2":os.path.join(rootDir,"software/miniconda3/envs/kyVirus/bin/bowtie2"),
"bin_mpa_v30_CHOCOPhlAn_201901":os.path.join(rootDir,"kyVirus","data","metaphlan_databases","mpa_v30_CHOCOPhlAn_201901"),
"bin_pythonVirus":os.path.join(rootDir,"software/miniconda3/envs/kyVirus/bin/python"),
"bin_metaphlan2":os.path.join(rootDir,"software/miniconda3/envs/kyVirus/bin/metaphlan2.py"),
"bin_metaphlanDatabases":os.path.join(rootDir,"kyVirus","data","metaphlan_databases"),
"bin_lib":os.path.join(rootDir,"software/miniconda3/envs/kyVirus/lib"),
"bin_kraken2":os.path.join(rootDir,"kyVirus","software","kraken2","kraken2"),
"bin_kraken2DB":os.path.join(rootDir,"kyVirus","data","kraken2_DB"),
"bin_bracken":os.path.join(rootDir,"kyVirus","software","bracken","bracken"),
"bin_value1":"\"SLIDINGWINDOW:4:20 MINLEN:35\"",
"bin_value2":"\"--very-sensitive --dovetail\"",
"bin_merge_metaphlan_tables":os.path.join(rootDir,"kyVirus/miniconda3/envs/kyVirus/bin/merge_metaphlan_tables.py"),#/mnt/nas/pipelines/CNVseq100k/pipline/kyVirus/miniconda3/envs/kyVirus/bin//merge_metaphlan_tables.py
"bin_combine_bracken_outputs":os.path.join(rootDir,"kyVirus/software/bracken/analysis_scripts/combine_bracken_outputs.py"),
"bin_find_positive":os.path.join(rootDir,"kyVirus/pathogeny/find_positive.py"),
"sc_virus":os.path.join(rootDir,"script","virus.wdl"),
"sc_virus_find_positive":os.path.join(rootDir,"script","virus.find_positive.wdl"),
"glibc_ky":"/home/bgiky/env/glibc/lib64/",
"glibc_sys":"/lib64/",
"bin_virusSampels":os.path.join(rootDir,"script","virus_samples.py"),
"bin_virusBracken":os.path.join(rootDir,"script","virus_bracken.py"),
"knead_sh":os.path.join(rootDir,"script","knead_outfile.sh")
}



Shell = os.path.join(work_path, "Shell")
Json = os.path.join(work_path, "Shell", "json")
if not os.path.exists(Shell):os.mkdir(Shell)
if not os.path.exists(Json):os.mkdir(Json)


##采用多次for循环 虽然不是最快的  但是比较清晰
######################soap.json and soap.sh##########################
"""
{
"wf_soap.fqPath":"/ifs9/zebra/MGISEQ-2000/V100400190292/V350086491/L04/V350086491_L04_41.fq.gz",
"wf_soap.sampleName":"22S03625076",
"wf_soap.sampleLib":"22L03625076-41",
"wf_soap.sampleBarcode":"41",
"wf_soap.seqPlatform":"MGISEQ-2000",
"wf_soap.outdir":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/test",
"wf_soap.bin_soap":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/software/SOAP",
"wf_soap.bin_hg19Index":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/SE/Hg19_reference/human.fa.index"
}
"""
"""
/home/KYpub/miniconda3/envs/kyVirus/bin/java -jar /dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/other/cromwell.jar run /dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/script/soap.wdl -i /dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/test/soap.json
"""

##sampleLib seqPlatform barcode fqpath fqpath2
for sampleID in fqpathDict.keys():
	sampleDict = dict()
	#sampleDict["wf_soap.fq1Path"] = fqpathDict[sampleID][3]
	sampleDict["wf_soap.sampleID"] = sampleID
	sampleDict["wf_soap.sampleLib"] = fqpathDict[sampleID][0]
	sampleDict["wf_soap.barcode"] = fqpathDict[sampleID][2]
	sampleDict["wf_soap.seqPlatform"] = fqpathDict[sampleID][1]
	sampleDict["wf_soap.outdir"] = work_path
	sampleDict["wf_soap.bin_soap"] = "{bin_soap}".format(**arg_dic)
	sampleDict["wf_soap.bin_hg19Index"] = "{bin_hg19Index}".format(**arg_dic)
	sampleSoapJson = os.path.join(Json,"Step1-1.soap."+sampleID+".json")
	soapShell = os.path.join(Shell,"Step1-1.soap."+sampleID+".sh")
	with open(sampleSoapJson,"w") as f:
		json.dump(sampleDict,f)
	f.close()
	with open(soapShell,"w") as ff:
		ff.write('#!/bin/bash\nif [ ! -d \"%s\" ]; then\n  mkdir -p %s\nfi\n'%(sampleDict["wf_soap.outdir"] ,sampleDict["wf_soap.outdir"] ))
		ff.write('cd %s/Shell\n'%(sampleDict["wf_soap.outdir"]))
		ff.write("{bin_java} -jar {bin_cromwell} run {sc_soap}".format(**arg_dic)+" -i "+sampleSoapJson+"\n")
	ff.close()
	"""
	if int(fqpathDict[sampleID][4]) != 0:
		sampleDict["wf_soap.fqPath"] = fqpathDict[sampleID][4]
		sampleDict["wf_soap.seqPlatform"] = fqpathDict[sampleID][1] + "_1"
		sampleDict["wf_soap.outdir"] = work_path
		sampleSoap2Json = os.path.join(Json,"Step1-1-0.soap."+sampleID+".json")
		with open(sampleSoap2Json,"w") as f:
			json.dump(sampleDict,f)
		f.close()
		cmd = '#!/bin/bash\nif [ ! -d \"%s\" ]; then\n  mkdir -p %s\nfi\n'%(sampleDict["wf_soap.outdir"],sampleDict["wf_soap.outdir"])
		cmd = cmd +"{bin_java} -jar {bin_cromwell} run {sc_soap}".format(**arg_dic)+" -i "+sampleSoap2Json+"\n"
		open(soapShell,"a+").write(cmd)
	"""

##############################soap end###############################

#####################fastq-1 start###################################
"""
{
"wf_qc.fqPath":"/ifs9/zebra/MGISEQ-2000/V100400190292/V350086491/L04/V350086491_L04_41.fq.gz",
"wf_qc.sampleName":"22S03625076",
"wf_qc.sampleLib":"22L03625076-41",
"wf_qc.sampleBarcode":"41",
"wf_qc.seqPlatform":"MGISEQ-2000",
"wf_qc.outdir":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/test/Hg19",
"wf_qc.bin_fqsta":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/SE/fq.sta.pl",
"wf_qc.bin_hg19Index":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/SE/Hg19_reference/human.fa.index"
}

"""
"""
/home/KYpub/miniconda3/envs/kyVirus/bin/java -jar /dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/other/cromwell.jar run /dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/script/fastq-1.wdl -i /dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/test/fastq-1.json

"""
###sampleLib seqPlatform barcode fqpath fqpath2
for sampleID in fqpathDict.keys():
	sampleDict = dict()
	#sampleDict["wf_qc.fq1Path"]=fqpathDict[sampleID][3]
	sampleDict["wf_qc.sampleID"]=sampleID
	sampleDict["wf_qc.sampleLib"]=fqpathDict[sampleID][0]
	sampleDict["wf_qc.barcode"]=fqpathDict[sampleID][2]
	sampleDict["wf_qc.seqPlatform"]=fqpathDict[sampleID][1]
	sampleDict["wf_qc.outdir"]=work_path
	sampleDict["wf_qc.bin_fqsta"]="{bin_fqsta}".format(**arg_dic)
	sampleDict["wf_qc.bin_perl"]="{bin_perl}".format(**arg_dic)
	sampleDict["wf_qc.bin_hg19Index"]="{bin_hg19Index}".format(**arg_dic)
	samplefastaJson = os.path.join(Json,"Step1-1.fastq1."+sampleID+".json")
	fastq1Shell = os.path.join(Shell,"Step1-1.fastq1."+sampleID+".sh")
	with open(samplefastaJson,"w") as f:
		json.dump(sampleDict,f)
	f.close()
	outdir = sampleDict["wf_qc.outdir"] + "/Hg19/" + sampleID + "/" + fqpathDict[sampleID][0] + "/" + fqpathDict[sampleID][2]+"_"+sampleDict["wf_qc.seqPlatform"]
	with open(fastq1Shell,"w") as ff:
		ff.write('#!/bin/bash\nif [ ! -d \"%s\" ]; then\n  mkdir -p %s\nfi\n'%(outdir,outdir))
		ff.write('cd %s/Shell\n'%(sampleDict["wf_qc.outdir"]))		
		ff.write("{bin_java} -jar {bin_cromwell} run {sc_fastq1}".format(**arg_dic)+" -i "+samplefastaJson+"\n")
	ff.close()
	if int(fqpathDict[sampleID][4]) != 0:
		sampleDict["wf_qc.fq1Path"]=fqpathDict[sampleID][4]
		sampleDict["wf_qc.seqPlatform"]=fqpathDict[sampleID][1] + "_1"
		samplefasta12Json = os.path.join(Json,"Step1-1-0.fastq1."+sampleID+".json")
		with open(samplefasta12Json,"w") as f:
			json.dump(sampleDict,f)
		f.close()
		outdir = sampleDict["wf_qc.outdir"] + "/Hg19/" + sampleID + "/" + fqpathDict[sampleID][0] + "/" + fqpathDict[sampleID][2]+"_"+sampleDict["wf_qc.seqPlatform"]+"_"
		cmd ='#!/bin/bash\nif [ ! -d \"%s\" ]; then\n  mkdir -p %s\nfi\n'%(outdir,outdir)
		cmd = cmd + "{bin_java} -jar {bin_cromwell} run {sc_fastq1}".format(**arg_dic)+" -i "+samplefasta12Json+"\n"
		open(fastq1Shell, "a+").write(cmd)



#############################fastq end###############################

#############################rmdup start#############################
"""
{
"wf_rmdup.gzsoap":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/test/Hg19/22S03625076/22L03625076-41/41_MGISEQ-2000/22S03625076_22L03625076-41_41_MGISEQ-2000.soap.gz",
"wf_rmdup.outdirStr":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/test/Hg19/22S03625076/22L03625076-41/41_MGISEQ-2000/22S03625076_22L03625076-41_41_MGISEQ-2000",
"wf_rmdup.bin_perl":"perl",
"wf_rmdup.bin_rmdup":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/SE/rmDup.pl",
"wf_rmdup.bin_SOAPcoverage":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/software/SOAPcoverage",
"wf_rmdup.bin_hg19":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/SE/Hg19_reference/human.fa",
"wf_rmdup.bin_Nregion":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/SE/HG19.N_region"
}
"""
"""
/home/KYpub/miniconda3/envs/kyVirus/bin/java -jar /dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/other/cromwell.jar run /dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/script/rmdup.wdl -i /dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/test/rmdup.json
"""
###sampleLib seqPlatform barcode fqpath fqpath2
for sampleID in fqpathDict.keys():
	sampleDict = dict()
	outdirStr = work_path+"/Hg19/"+sampleID+"/"+fqpathDict[sampleID][0]+"/"+fqpathDict[sampleID][2]+"_"+fqpathDict[sampleID][1]+"/"+sampleID+"_"+fqpathDict[sampleID][0]+"_"+fqpathDict[sampleID][2]+"_"+fqpathDict[sampleID][1]
	sampleDict["wf_rmdup.sampleID"] = sampleID
	sampleDict["wf_rmdup.outdir"] = work_path
	#sampleDict["wf_rmdup.gzsoap"] = outdirStr + ".soap.gz"
	sampleDict["wf_rmdup.outdir_hg19"] = outdirStr
	sampleDict["wf_rmdup.bin_perl"] = "{bin_perl}".format(**arg_dic)
	sampleDict["wf_rmdup.bin_rmdup"] = "{bin_rmdup}".format(**arg_dic)
	sampleDict["wf_rmdup.bin_SOAPcoverage"] = "{bin_SOAPcoverage}".format(**arg_dic)
	sampleDict["wf_rmdup.bin_hg19"] = "{bin_hg19}".format(**arg_dic)
	sampleDict["wf_rmdup.bin_Nregion"] = "{bin_Nregion}".format(**arg_dic)
	samplermdupJson = os.path.join(Json,"Step1-2.rmdup."+sampleID+".json")
	rmdupShell = os.path.join(Shell,"Step1-2.rmdup."+sampleID+".sh")
	with open(samplermdupJson,"w") as f:
		json.dump(sampleDict,f)
	f.close()
	with open(rmdupShell,"w") as ff:
		ff.write('#!/bin/bash\n')
		ff.write('cd %s/Shell\n'%(sampleDict["wf_rmdup.outdir"]))
		ff.write("{bin_java} -jar {bin_cromwell} run {sc_rmdup}".format(**arg_dic)+" -i "+samplermdupJson+"\n")
	ff.close()
	"""
	if int(fqpathDict[sampleID][4]) != 0: 
		outdirStr = work_path+"/Hg19/"+sampleID+"/"+fqpathDict[sampleID][0]+"/"+fqpathDict[sampleID][2]+"_"+fqpathDict[sampleID][1]+"_1"+"/"+sampleID+"_"+fqpathDict[sampleID][0]+"_"+fqpathDict[sampleID][2]+"_"+fqpathDict[sampleID][1]+"_1"
		sampleDict["wf_rmdup.gzsoap"] = outdirStr + ".soap.gz"
		sampleDict["wf_rmdup.outdirStr"] = outdirStr
		samplermdup2Json = os.path.join(Json,"Step1-2-0.rmdup."+sampleID+".json")
		with open(samplermdup2Json,"w") as f:
			json.dump(sampleDict,f)
		f.close()
		cmd = "{bin_java} -jar {bin_cromwell} run {sc_rmdup}".format(**arg_dic)+" -i "+samplermdup2Json+"\n"
		open(rmdupShell, 'a+').write(cmd)
	"""
#############################rmdup end###############################

#####################fastq-1 start###################################
"""
{
"wf_qc.gzrmdup":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/test/Hg19/22S03625076/22L03625076-41/41_MGISEQ-2000/22S03625076_22L03625076-41_41_MGISEQ-2000.soap.rmdup.gz",
"wf_qc.outInfo":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/test/Hg19/22S03625076/22L03625076-41/41_MGISEQ-2000/22S03625076_22L03625076-41_41_MGISEQ-2000.info",
"wf_qc.bin_perl":"perl",
"wf_qc.bin_soapSta":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/SE/cov.soap.sta.pl",
"wf_qc.dupStat":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/test/Hg19/22S03625076/22L03625076-41/41_MGISEQ-2000/22S03625076_22L03625076-41_41_MGISEQ-2000.rmdup.stat",
"wf_qc.cvgCoverageInfo":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/test/Hg19/22S03625076/22L03625076-41/41_MGISEQ-2000/22S03625076_22L03625076-41_41_MGISEQ-2000.cvg_coverage.info"}
"""
"""
/home/KYpub/miniconda3/envs/kyVirus/bin/java -jar /dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/other/cromwell.jar run /dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/script/fastq-2.wdl -i /dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/test/fastq-2.json
"""
for sampleID in fqpathDict.keys():
	sampleDict = dict()
	outdirStr = work_path+"/Hg19/"+sampleID+"/"+fqpathDict[sampleID][0]+"/"+fqpathDict[sampleID][2]+"_"+fqpathDict[sampleID][1]+"/"+sampleID+"_"+fqpathDict[sampleID][0]+"_"+fqpathDict[sampleID][2]+"_"+fqpathDict[sampleID][1]
	#sampleDict["wf_qc.gzrmdup"] = outdirStr + ".soap.rmdup.gz"
	sampleDict["wf_qc.outdir_hg19"] = outdirStr
	sampleDict["wf_qc.out_Info"] = outdirStr + ".info"
	sampleDict["wf_qc.bin_perl"] = "{bin_perl}".format(**arg_dic)
	sampleDict["wf_qc.bin_soapSta"] = "{bin_soapSta}".format(**arg_dic)
	#sampleDict["wf_qc.dupStat"] = outdirStr + ".rmdup.stat"
	#sampleDict["wf_qc.cvgCoverageInfo"] = outdirStr + ".cvg_coverage.info"
	samplefastq2Json = os.path.join(Json,"Step1-3.fastq2."+sampleID+".json")
	fastq2Shell = os.path.join(Shell,"Step1-3.fastq2."+sampleID+".sh")
	with open(samplefastq2Json,"w") as f:
		json.dump(sampleDict,f)
	f.close()
	with open(fastq2Shell,"w") as ff:
		ff.write('#!/bin/bash\n')
		ff.write('cd %s/Shell\n'%(work_path))
		ff.write("{bin_java} -jar {bin_cromwell} run {sc_fastq2}".format(**arg_dic)+" -i "+samplefastq2Json+"\n")
	ff.close()
	if int(fqpathDict[sampleID][4]) != 0:
		outdirStr = work_path+"/Hg19/"+sampleID+"/"+fqpathDict[sampleID][0]+"/"+fqpathDict[sampleID][2]+"_"+fqpathDict[sampleID][1]+"_1"+"/"+sampleID+"_"+fqpathDict[sampleID][0]+"_"+fqpathDict[sampleID][2]+"_"+fqpathDict[sampleID][1]+"_1"
		sampleDict["wf_qc.gzrmdup"] = outdirStr + ".soap.rmdup.gz"
		sampleDict["wf_qc.outInfo"] = outdirStr + ".info"
		sampleDict["wf_qc.dupStat"] = outdirStr + ".rmdup.stat"
		samplefastq20Json = os.path.join(Json,"Step1-3.fastq2."+sampleID+".json")
		with open(samplefastq20Json,"w") as f:
			json.dump(sampleDict,f)
		f.close()
		cmd = "{bin_java} -jar {bin_cromwell} run {sc_fastq2}".format(**arg_dic)+" -i "+samplefastq20Json+"\n"
		oepn(fastq2Shell,"a+").write(cmd)


		
#####################fastq-1 end#####################################

#####################DataProduction start############################
"""
{
"wf_DataProduction.InformationList":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/test/Hg19/Information.list",
"wf_DataProduction.outFile":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/test/Hg19/DataProduction.xls",
"wf_DataProduction.bin_perl":"perl",
"wf_DataProduction.bin_DataProduction":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/SE/Data_production.pl"
}
"""
"""
/home/KYpub/miniconda3/envs/kyVirus/bin/java -jar /dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/other/cromwell.jar run /dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/script/Data_Production.wdl -i /dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/test/Data_Production.json
"""
sampleDict = dict()
sampleDict["wf_DataProduction.bin_perl"] = "{bin_perl}".format(**arg_dic)
sampleDict["wf_DataProduction.bin_DataProduction"] = "{bin_DataProduction}".format(**arg_dic)
sampleDict["wf_DataProduction.outdir"] = work_path
#sampleDict["wf_DataProduction.bin_python"] = "{bin_python}".format(**arg_dic)
#sampleDict["wf_DataProduction.bin_pscc"] = "{bin_pscc}".format(**arg_dic)
DataProductionJson = os.path.join(Json,"Step1-4.DataProduction.json")
DataProductionShell = os.path.join(Shell,"Step1-4.DataProduction.sh")
with open(DataProductionJson,"w") as f:
	json.dump(sampleDict,f)
f.close()
with open(DataProductionShell,"w") as ff:
	ff.write('#!/bin/bash\n')
	ff.write('cd %s/Shell\n'%(work_path))
	ff.write("{bin_java} -jar {bin_cromwell} run {sc_DataProduction}".format(**arg_dic)+" -i "+DataProductionJson+"\n"    )
ff.close()

####################DataProduction end###############################

#####################pscc start########################################
for sampleID in fqpathDict.keys():
	sampleDict = dict()
	sampleDict["wf_pscc.gzrmdup"] = work_path+"/Hg19/"+sampleID+"/"+fqpathDict[sampleID][0]+"/"+fqpathDict[sampleID][2]+"_"+fqpathDict[sampleID][1]+"/"+sampleID+"_"+fqpathDict[sampleID][0]+"_"+fqpathDict[sampleID][2]+"_"+fqpathDict[sampleID][1]+".soap.rmdup.gz"
	if int(fqpathDict[sampleID][4]) != 0:
		sampleDict["wf_pscc.gzrmdup"] = sampleDict["wf_pscc.gzrmdup"] + work_path+"/Hg19/"+sampleID+"/"+fqpathDict[sampleID][0]+"/"+fqpathDict[sampleID][2]+"_"+fqpathDict[sampleID][1]+"_1/"+sampleID+"_"+fqpathDict[sampleID][0]+"_"+fqpathDict[sampleID][2]+"_"+fqpathDict[sampleID][1]+"_1.soap.rmdup.gz"
	sampleDict["wf_pscc.outdir"] = work_path
	sampleDict["wf_pscc.sampleID"] = sampleID
	sampleDict["wf_pscc.bin_perl"] = "{bin_perl}".format(**arg_dic)
	sampleDict["wf_pscc.bin_hg19"] = "{bin_hg19}".format(**arg_dic)
	sampleDict["wf_pscc.bin_soap2ext"] = "{bin_soap2ext}".format(**arg_dic)
	sampleDict["wf_pscc.bin_tagsgc"] = "{bin_tagsgc}".format(**arg_dic)
	sampleDict["wf_pscc.bin_Aneuploid"] = "{bin_Aneuploid}".format(**arg_dic)
	sampleDict["wf_pscc.bin_winRatio"] = "{bin_winRatio}".format(**arg_dic) + infoDict[sampleID][2]
	sampleDict["wf_pscc.bin_RatioLoose"] = "{bin_RatioLoose}".format(**arg_dic)
	sampleDict["wf_pscc.bin_Segmentation"] = "{bin_Segmentation}".format(**arg_dic)
	sampleDict["wf_pscc.bin_cnvFilter"] = "{bin_cnvFilter}".format(**arg_dic)
	sampleDict["wf_pscc.bin_hg19Length"] = "{bin_hg19Length}".format(**arg_dic)
	sampleDict["wf_pscc.bin_hg19Nregion"] = "{bin_hg19Nregion}".format(**arg_dic)
	sampleDict["wf_pscc.bin_spValue"] = "{bin_spValue}".format(**arg_dic)
	sampleDict["wf_pscc.bin_ppValue"] = "{bin_ppValue}".format(**arg_dic)
	sampleDict["wf_pscc.bin_Graph"] = "{bin_Graph}".format(**arg_dic)
	sampleDict["wf_pscc.gender_sh"] = "{gender_sh}".format(**arg_dic)
	sampleDict["wf_pscc.glibc_ky"] = "{glibc_ky}".format(**arg_dic)
	sampleDict["wf_pscc.glibc_sys"] = "{glibc_sys}".format(**arg_dic)
	sampleDict["wf_pscc.bin_GraphReport"] = "{bin_GraphReport}".format(**arg_dic)
	samplePSCCJson = os.path.join(Json,"Step2-1.sample."+sampleID+".json")
	psccShell = os.path.join(Shell,"Step2-1.sample."+sampleID+".sh")
	outdir = work_path + "/Hg19.cnv/Case/" + sampleID
	outdirplot = work_path + "/Hg19.cnv/Case/" + sampleID + "/" + sampleID + ".plot"
	with open(samplePSCCJson,"w") as f:
		json.dump(sampleDict,f)
	f.close()
	with open(psccShell,"w") as ff:
		ff.write('#!/bin/bash\nif [ ! -d \"%s\" ]; then\n  mkdir -p %s\nfi\n'%(outdir,outdir))
		ff.write('if [ ! -d \"%s\" ]; then\n  mkdir -p %s\nfi\n'%(outdirplot,outdirplot))
		ff.write('cd %s/Shell\n'%(work_path))
		ff.write("{bin_java} -jar {bin_cromwell} run {sc_pscc}".format(**arg_dic)+" -i "+samplePSCCJson+"\n")
	ff.close()
####################pscc end########################################


#####################GC_correction_ratio start#######################
"""
{
"wf_GCcorrection.gzext":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/test/Hg19.cnv/Case/22S03625076/22S03625076.ext.gz",
"wf_GCcorrection.gzratio":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/test/Hg19.cnv/Case/22S03625076/22S03625076.ratio.gz",
"wf_GCcorrection.outdir":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/test/Hg19.cnv/Slope",
"wf_GCcorrection.sampleName":"22S03625076",
"wf_GCcorrection.bin_perl":"perl",
"wf_GCcorrection.bin_EXT_SPLIT":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/SLOPE-code-use-now/bin/EXT_SPLIT_individual_v1.pl",
"wf_GCcorrection.bin_slopeWinRatio":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/SLOPE-code-use-now/control_cnv/winRatio_MGISEQ-2000",
"wf_GCcorrection.bin_GCcorrection":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/SLOPE-code-use-now/bin/GC_correction.pl",
"wf_GCcorrection.bin_splitRatio":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/SLOPE-code-use-now/bin/split_for_ratio.pl"
}
"""
"""
/home/KYpub/miniconda3/envs/kyVirus/bin/java -jar /dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/other/cromwell.jar run /dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/script/GC_correction_ratio.wdl -i /dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/test/GC_correction_ratio.json
"""
###sampleLib seqPlatform barcode fqpath fqpath2
for sampleID in fqpathDict.keys():
	sampleDict = dict()
	sampleDict["wf_GCcorrection.outdir"] = work_path 
	sampleDict["wf_GCcorrection.sampleID"] = sampleID
	sampleDict["wf_GCcorrection.bin_perl"] = "{bin_perl}".format(**arg_dic)
	sampleDict["wf_GCcorrection.bin_EXT_SPLIT"] = "{bin_EXT_SPLIT}".format(**arg_dic)
	sampleDict["wf_GCcorrection.bin_slopeWinRatio"] = "{bin_slopeWinRatio}".format(**arg_dic) + infoDict[sampleID][2]
	sampleDict["wf_GCcorrection.bin_GCcorrection"] = "{bin_GCcorrection}".format(**arg_dic)
	sampleDict["wf_GCcorrection.bin_splitRatio"] = "{bin_splitRatio}".format(**arg_dic)
	sampleGCRatioJson = os.path.join(Json,"Step2-2.GC_correction_ratio."+sampleID+".json")
	GCRatioShell = os.path.join(Shell,"Step2-2.GC_correction_ratio."+sampleID+".sh")
	with open(sampleGCRatioJson,"w") as f:
		json.dump(sampleDict,f)
	f.close()
	gcoutdir = sampleDict["wf_GCcorrection.outdir"] + "/" + sampleID + "/" + sampleID + "_GC"
	ratiooutdir = sampleDict["wf_GCcorrection.outdir"] + "/" + sampleID + "/" + sampleID + "_ratio"
	with open(GCRatioShell,"w") as ff:
		ff.write('#!/bin/bash\nif [ ! -d \"%s\" ]; then\n  mkdir -p %s\nfi\n'%(gcoutdir,gcoutdir))
		ff.write('if [ ! -d \"%s\" ]; then\n  mkdir %s\nfi\n'%(ratiooutdir,ratiooutdir))
		ff.write('cd %s/Shell\n'%(work_path))
		ff.write("{bin_java} -jar {bin_cromwell} run {sc_GC_correction_ratio}".format(**arg_dic)+" -i "+sampleGCRatioJson+"\n")
	ff.close()

#####################GC_correction_ratio end#########################

#####################Slope start#######################################
"""
{
"wf_slope.outdir":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/test",
"wf_slope.sampleID":"22S03625076",
"wf_slope.bin_perl":"perl",
"wf_slope.bin_zibizslope":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/SLOPE-code-use-now/bin/zibiz-slope_P_v18.pl",
"wf_slope.bin_find_merge":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/SLOPE-code-use-now/bin/find_merge_v2.pl",
"wf_slope.bin_CNVFilter":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/SLOPE-code-use-now/bin/CNVFilter.pl"
}
"""
for sampleID in fqpathDict.keys():
	sampleDict = dict()
	sampleDict["wf_slope.outdir"] = work_path
	sampleDict["wf_slope.sampleID"] = sampleID
	sampleDict["wf_slope.bin_perl"] = "{bin_perl}".format(**arg_dic)
	sampleDict["wf_slope.bin_zibizslope"] = "{bin_zibizslope}".format(**arg_dic)
	sampleDict["wf_slope.bin_find_merge"] = "{bin_find_merge}".format(**arg_dic)
	sampleDict["wf_slope.bin_CNVFilter"] = "{bin_CNVFilter}".format(**arg_dic)
	sampleDict["wf_slope.gender_sh"] = "{gender_sh}".format(**arg_dic)
	sampleSlopeJson = os.path.join(Json,"Step2-3.slope."+sampleID+".json")
	slopeShell = os.path.join(Shell,"Step2-3.slope."+sampleID+".sh")
	with open(sampleSlopeJson,"w") as f:
		json.dump(sampleDict,f)
	f.close()
	with open(slopeShell,"w") as ff:
		ff.write('#!/bin/bash\n')
		ff.write('cd %s/Shell\n'%(work_path))
		ff.write("{bin_java} -jar {bin_cromwell} run {sc_slope}".format(**arg_dic)+" -i "+sampleSlopeJson+    "\n")
	ff.close()
###################Slope end#########################################

#####################step4 start######################################
"""
{
"wf_highConfidenceCNV.outdir":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/test",
"wf_highConfidenceCNV.sampleName":"22S03625076",
"wf_highConfidenceCNV.productID":"DX0007",
"wf_highConfidenceCNV.bin_python":"/dfsyt1/B2C_COM_S1/B2C_PSCC/pipline/kyCNV20211022/pipline/AIO_sourcecode/tools/python",
"wf_highConfidenceCNV.bin_gencheck":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/AIO_sourcecode/cnv/gencheck.py",
"wf_highConfidenceCNV.bin_merge_5M_pro":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/AIO_sourcecode/result/merge_5M_pro.py",
"wf_highConfidenceCNV.bin_merge100K":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/AIO_sourcecode/result/merge_100K_v1.py",
"wf_highConfidenceCNV.bin_highConfidence100K":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/AIO_sourcecode/result/high_confidence_100K.py",
"wf_highConfidenceCNV.bin_hg19Length":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/PSCC/Attachment/hg19.length",
"wf_highConfidenceCNV.bin_hg19cytoBand":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/PSCC/Attachment/hg19.cytoBand",
"wf_highConfidenceCNV.bin_AIOdata":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/AIO_sourcecode/dat"
}
"""
#/home/KYpub/miniconda3/envs/kyVirus/bin/java -jar /dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/other/cromwell.jar run /dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/script/step4.wdl -i /dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/test/step4.json
###sampleLib seqPlatform barcode fqpath fqpath2
##gender 100K/1M hospitalName 
for sampleID in fqpathDict.keys():
	sampleDict = dict()
	sampleDict["wf_highConfidenceCNV.outdir"] = work_path
	sampleDict["wf_highConfidenceCNV.sampleID"] = sampleID
	sampleDict["wf_highConfidenceCNV.productID"] = infoDict[sampleID][3]
	sampleDict["wf_highConfidenceCNV.bin_python"] = "{bin_python}".format(**arg_dic)
	sampleDict["wf_highConfidenceCNV.bin_gencheck"] = "{bin_gencheck}".format(**arg_dic)
	sampleDict["wf_highConfidenceCNV.bin_merge_5M_pro"] = "{bin_merge_5M_pro}".format(**arg_dic)
	sampleDict["wf_highConfidenceCNV.bin_merge100K"] = "{bin_merge100K}".format(**arg_dic)
	sampleDict["wf_highConfidenceCNV.bin_hg19Nregion"] = "{bin_hg19Nregion}".format(**arg_dic) 
	sampleDict["wf_highConfidenceCNV.bin_highConfidence100K"] = "{bin_highConfidence100K}".format(**arg_dic)
	sampleDict["wf_highConfidenceCNV.bin_hg19Length"] = "{bin_hg19Length}".format(**arg_dic)
	sampleDict["wf_highConfidenceCNV.bin_hg19cytoBand"] = "{bin_hg19cytoBand}".format(**arg_dic)
	sampleDict["wf_highConfidenceCNV.bin_AIOdata"] = "{bin_AIOdata}".format(**arg_dic)
	sampleDict["wf_highConfidenceCNV.product_sh"] = "{product_sh}".format(**arg_dic)
	sampleStep4Json = os.path.join(Json,"Step2-4."+sampleID+".json")
	step4Shell = os.path.join(Shell,"Step2-4."+sampleID+".sh")
	with open(sampleStep4Json,"w") as f:
		json.dump(sampleDict,f)
	f.close()
	with open(step4Shell,"w") as ff:
		ff.write('#!/bin/bash\n')
		ff.write('cd %s/Shell\n'%(work_path))
		ff.write("{bin_java} -jar {bin_cromwell} run {sc_step4}".format(**arg_dic)+" -i "+sampleStep4Json+"\n")
	ff.close()


#####################step4 end########################################

######################result start####################################
"""
{
"wf_result.outdir":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/test",
"wf_result.bin_python3":"/home/KYpub/miniconda3/envs/kyCNV/bin/python",
"wf_result.bin_python2":"/home/KYpub/miniconda3/envs/py2/bin/python",
"wf_result.bin_compareResult":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/AIO_sourcecode/result/compare_result_pro.py",
"wf_result.bin_mergepng":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/AIO_sourcecode/result/pscc_merge_png.py"
}
"""
#/home/KYpub/miniconda3/envs/kyVirus/bin/java -jar /dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/other/cromwell.jar run /dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/script/result.wdl -i /dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/test/result.json
###sampleLib seqPlatform barcode fqpath fqpath2
##gender 100K/1M hospitalName
sampleDict = dict()
sampleDict["wf_result.outdir"] = work_path
sampleDict["wf_result.bin_python3"] = "{bin_python}".format(**arg_dic)
sampleDict["wf_result.bin_python2"] = "{bin_python2}".format(**arg_dic)
sampleDict["wf_result.bin_compareResult"] = "{bin_compareResult}".format(**arg_dic)
sampleDict["wf_result.bin_mergepng"] = "{bin_mergepng}".format(**arg_dic)
sampleResultJson = os.path.join(Json,"Step3.result.json")
resultShell = os.path.join(Shell,"Step3.result.sh")
with open(sampleResultJson,"w") as f:
	json.dump(sampleDict,f)
f.close()
outResult = work_path + "/result"
outResultpng = outResult + "/for_check"
with open(resultShell,"w") as ff:
	ff.write('#!/bin/bash\nif [ ! -d \"%s\" ]; then\n  mkdir %s\nfi\nif [ -d  \"%s\" ]; then \n rm -r %s\nfi\n'%(outResult,outResult,outResultpng,outResultpng))
	ff.write('cd %s/Shell\n'%(work_path))
	ff.write("{bin_java} -jar {bin_cromwell} run {sc_result}".format(**arg_dic)+" -i "+sampleResultJson+"\n")
ff.close()


######################result end######################################
"""
######################Annotation start################################
{
"wf_Annotation.resultFile":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/test/result/result.xlsx",
"wf_Annotation.bin_python":"/home/KYpub/miniconda3/envs/kyCNV/bin/python",
"wf_Annotation.bin_Annotation":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/Annotation_V3.0/KY_Annotation.py"
}
#/home/KYpub/miniconda3/envs/kyVirus/bin/java -jar /dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/other/cromwell.jar run /dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/script/Anno.wdl -i /dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/test/Anno.json

sampleDict = dict()
sampleDict["wf_Annotation.resultFile"] = work_path + "/result/result.xlsx" 
sampleDict["wf_Annotation.bin_python"] = "{bin_python}".format(**arg_dic)
sampleDict["wf_Annotation.bin_Annotation"] = "{bin_Annotation}".format(**arg_dic)
sampleAnnoJson = os.path.join(Json,"Step4.Anno.json")
annoShell = os.path.join(Shell,"Step4.Anno.sh")
with open(sampleAnnoJson,"w") as f:
	json.dump(sampleDict,f)
f.close()
with open(annoShell,"w") as ff:
	ff.write("{bin_java} -jar {bin_cromwell} run {sc_anno}".format(**arg_dic)+" -i "+sampleAnnoJson+"\n")
ff.close()
######################Annotation end##################################
"""

######################virus start#####################################
"""
{
"wf_virus.fqpath":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/test/V350086491_L04_41.fq.gz",
"wf_virus.outdir":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/test",
"wf_virus.bin_kneaddata":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/kyVirus/software/kyVirus/bin/kneaddata",
"wf_virus.bin_kneaddataDatabase":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/kyVirus/data/kneaddata_database",
"wf_virus.bin_trimmomatic":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/kyVirus/software/kyVirus/share/trimmomatic-0.39-1/",
"wf_virus.bin_software":"/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/kyVirus/software/kyVirus/bin",
"wf_virus.sampleID": "22S03625076",
"wf_virus.bin_bowtie2": "/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/kyVirus/software/kyVirus/bin/bowtie2",
"wf_virus.bin_mpa_v30_CHOCOPhlAn_201901": "/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/kyVirus/data/metaphlan_databases/mpa_v30_CHOCOPhlAn_201901",
"wf_virus.bin_python": "/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/kyVirus/software/kyVirus/bin/python",
"wf_virus.bin_metaphlan2": "/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/kyVirus/software/kyVirus/bin/metaphlan2.py",
"wf_virus.bin_metaphlanDatabases": "/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/kyVirus/data/metaphlan_databases",
"wf_virus.bin_lib": "/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/kyVirus/software/kyVirus/lib",
"wf_virus.bin_kraken2": "/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/kyVirus/software/kraken2/kraken2",
"wf_virus.bin_kraken2DB": "/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/kyVirus/data/kraken2_DB",
"wf_virus.bin_bracken": "/dfsyt1/B2C_COM_S1/B2C_PSCC/user/manjianfen/WDL/bin/kyVirus/software/bracken/bracken",
"wf_virus.bin_value1":"\"SLIDINGWINDOW:4:20 MINLEN:35\"",
"wf_virus.bin_value2":"\"--very-sensitive --dovetail\""
}

"""
####sampleLib seqPlatform barcode fqpath fqpath2
##gender 100K/1M hospitalName

for sampleID in fqpathDict.keys():
	samplevirusJson = os.path.join(Json,"Step4-1.virus."+sampleID+".json")
	virusShell = os.path.join(Shell,"Step4-1.virus."+sampleID+".sh")
	if infoDict[sampleID][1] == "1M":
		os.system('touch {0}'.format(samplevirusJson))
		cmd = '#!/usr/bin/bash\n'
		open(os.path.join(Shell, 'Step4-1.virus.{0}.sh'.format(sampleID)),'w').write(cmd)
		continue
	sampleDict = dict()
	sampleDict["wf_virus.fq1path"] = fqpathDict[sampleID][3]
	sampleDict["wf_virus.outdir"] = work_path
	sampleDict["wf_virus.sampleID"] = sampleID
	sampleDict["wf_virus.barcode"] =  fqpathDict[sampleID][2]
	sampleDict["wf_virus.bin_software"] = "{bin_software}".format(**arg_dic)
	sampleDict["wf_virus.bin_kneaddata"] = "{bin_kneaddata}".format(**arg_dic)
	sampleDict["wf_virus.bin_kneaddataDatabase"] = "{bin_kneaddataDatabase}".format(**arg_dic)
	sampleDict["wf_virus.bin_trimmomatic"] = "{bin_trimmomatic}".format(**arg_dic)
	sampleDict["wf_virus.bin_bowtie2"] = "{bin_bowtie2}".format(**arg_dic)
	sampleDict["wf_virus.bin_mpa_v30_CHOCOPhlAn_201901"] = "{bin_mpa_v30_CHOCOPhlAn_201901}".format(**arg_dic)
	sampleDict["wf_virus.bin_python"] = "{bin_pythonVirus}".format(**arg_dic)
	sampleDict["wf_virus.bin_metaphlan2"] = "{bin_metaphlan2}".format(**arg_dic)
	sampleDict["wf_virus.bin_metaphlanDatabases"] = "{bin_metaphlanDatabases}".format(**arg_dic)
	sampleDict["wf_virus.bin_lib"] = "{bin_lib}".format(**arg_dic)
	sampleDict["wf_virus.bin_kraken2"] = "{bin_kraken2}".format(**arg_dic)
	sampleDict["wf_virus.bin_kraken2DB"] = "{bin_kraken2DB}".format(**arg_dic)
	sampleDict["wf_virus.bin_bracken"] = "{bin_bracken}".format(**arg_dic)
	sampleDict["wf_virus.bin_value1"] = "{bin_value1}".format(**arg_dic)
	sampleDict["wf_virus.bin_value2"] = "{bin_value2}".format(**arg_dic)
	with open(samplevirusJson,"w") as f:
		json.dump(sampleDict,f)
	f.close()
	with open(virusShell,"w") as ff:
		ff.write('#!/bin/bash\n')
		ff.write('cd %s/Shell\n'%(work_path))
		ff.write("{bin_java} -jar {bin_cromwell} run {sc_virus}".format(**arg_dic)+" -i "+samplevirusJson+"\n")
	ff.close()



######################virus end#######################################

######################virus report start##############################
sampleListVirus = list()
for sampleID in fqpathDict.keys():
	if infoDict[sampleID][1] == "100K":
		sampleListVirus.append(sampleID)
if len(sampleListVirus) == 0:
	virusReportShell = os.path.join(Shell,"Step4-2.virus.find_positive.sh")
	with open(virusReportShell,"w") as ff:
		ff.write('#!/bin/bash\n')
		#ff.write('cd %s/Shell\n'%(work_path))
		ff.write("touch "+work_path+"/result/virus_report.xlsx"+"\n")
		ff.close()
else:
	#samples = ",".join(sampleListVirus)
	#work_path_List = [work_path+"/virus/kraken_out/"] * len(sampleListVirus)
	#brackenxls = [work_path_List[i] + sampleListVirus[i] + ".bracken.xls"  for i in range(0, len(sampleListVirus))] 
	#xlsxs = " ".join(brackenxls)#test/virus/krakenOut/22B27012751.bracken.xls
	#/dfsyt1/B2C_COM_S1/B2C_PSCC/user/pangao/YTJ/Clinic_workspace/KY-20220623-TEST/10/test10-1/virus/knead_out/22S02604718_kneaddata.fastq.gz
	#work_path_List = [work_path+"/virus/knead_out/"] * len(sampleListVirus)
	#knead_outfile = [sampleListVirus[i] + " " + work_path_List[i] + sampleListVirus[i] + "_kneaddata.fastq.gz" for i in range(0, len(sampleListVirus))]
	sampleVirusReportJson = os.path.join(Json,"Step4-2.virus.find_positive.json")
	virusReportShell = os.path.join(Shell,"Step4-2.virus.find_positive.sh")
	sampleDict = dict()
	sampleDict["wf_virus_find_positive.outdir"] = work_path
	sampleDict["wf_virus_find_positive.bin_python"] = "{bin_python}".format(**arg_dic)
	sampleDict["wf_virus_find_positive.bin_merge_metaphlan_tables"] = "{bin_merge_metaphlan_tables}".format(**arg_dic)
	sampleDict["wf_virus_find_positive.bin_combine_bracken_outputs"] = "{bin_combine_bracken_outputs}".format(**arg_dic)
	sampleDict["wf_virus_find_positive.bin_find_positive"] = "{bin_find_positive}".format(**arg_dic)
	sampleDict["wf_virus_find_positive.knead_sh"] = "{knead_sh}".format(**arg_dic)
	sampleDict["wf_virus_find_positive.bin_virusSampels"] = "{bin_virusSampels}".format(**arg_dic)
	sampleDict["wf_virus_find_positive.bin_virusBracken"] = "{bin_virusBracken}".format(**arg_dic)

	with open(sampleVirusReportJson,"w") as f:
		json.dump(sampleDict,f)
	f.close()
	with open(virusReportShell,"w") as ff:
		ff.write('#!/bin/bash\n')
		ff.write('cd %s/Shell\n'%(work_path))
		ff.write("{bin_java} -jar {bin_cromwell} run {sc_virus_find_positive}".format(**arg_dic)+" -i "+sampleVirusReportJson+"\n")
	ff.close()
######################virus report end################################


##########################for qsub####################################
# 生成任务描述文件 #
logging.info('write allSteps.json')
allJsonList = list()
# Step1-1.Soap.19S6926253.sh
step1Dict = dict()
step1Dict['first'] = 1
step1Dict['name'] = 'Step1-1.soap'
step1Dict['priorStep'] = []
step1Dict['nextStep'] = ['Step1-1.fastq1']
jobShList = list()
for sampleId in fqpathDict.keys():
	jobDict = dict()
	jobDict['mem'] = '6'
	jobDict['sh'] = os.path.join(work_path, 'Shell', 'Step1-1.soap.{0}.sh'.format(sampleId))
	jobShList.append(jobDict)
step1Dict['jobSh'] = jobShList
allJsonList.append(step1Dict)
#Step1-1.fastq1
step1_1Dict = dict()
step1_1Dict['name'] = 'Step1-1.fastq1'
step1_1Dict['priorStep'] = ['Step1-1.soap']
step1_1Dict['nextStep'] = ['Step1-2.rmdup']
jobShList = list()
for sampleId in fqpathDict.keys():
	jobDict = dict()
	jobDict['mem'] = '1'
	jobDict['sh'] = os.path.join(work_path, 'Shell', 'Step1-1.fastq1.{0}.sh'.format(sampleId))
	jobShList.append(jobDict)
step1_1Dict['jobSh'] = jobShList
allJsonList.append(step1_1Dict)

#Step1-2.rmdup
step1_2Dict = dict()
step1_2Dict['name'] = 'Step1-2.rmdup'
step1_2Dict['priorStep'] = ['Step1-1.fastq1']
step1_2Dict['nextStep'] = ['Step1-3.fastq2']
jobShList = list()
for sampleId in fqpathDict.keys():
	jobDict = dict()
	jobDict['mem'] = '1'
	jobDict['sh'] = os.path.join(work_path, 'Shell', 'Step1-2.rmdup.{0}.sh'.format(sampleId))
	jobShList.append(jobDict)
step1_2Dict['jobSh'] = jobShList
allJsonList.append(step1_2Dict)
#Step1-3.fastq2

step1_3Dict = dict()
step1_3Dict['name'] = 'Step1-3.fastq2'
step1_3Dict['priorStep'] = ['Step1-2.rmdup']
step1_3Dict['nextStep'] = ['Step1-4.DataProduction']
jobShList = list()
for sampleId in fqpathDict.keys():
	jobDict = dict()
	jobDict['mem'] = '1'
	jobDict['sh'] = os.path.join(work_path, 'Shell', 'Step1-3.fastq2.{0}.sh'.format(sampleId))
	jobShList.append(jobDict)
step1_3Dict['jobSh'] = jobShList
allJsonList.append(step1_3Dict)



# Step1-4.DataProduction.sh
step1_4Dict = dict()
step1_4Dict['name'] = 'Step1-4.DataProduction'
step1_4Dict['priorStep'] = ['Step1-3.fastq3']
step1_4Dict['nextStep'] = ['Step2']
jobShList = list()
jobDict = dict()
jobDict['mem'] = '2'
jobDict['sh'] = os.path.join(work_path, 'Shell', 'Step1-4.DataProduction.sh')
jobShList.append(jobDict)
step1_4Dict['jobSh'] = jobShList
allJsonList.append(step1_4Dict)

# Step2-1.sample.19S6926253.sh
step2_1Dict = dict()
step2_1Dict['name'] = 'Step2-1.sample'
step2_1Dict['priorStep'] = ['Step1-2.DataProduction']
step2_1Dict['nextStep'] = ['Step2-2']
jobShList = list()
for sampleId in fqpathDict.keys():
	jobDict = dict()
	jobDict['mem'] = '2'
	jobDict['sh'] = os.path.join(work_path, 'Shell', 'Step2-1.sample.{0}.sh'.format(sampleId))
	jobShList.append(jobDict)
step2_1Dict['jobSh'] = jobShList
allJsonList.append(step2_1Dict)

#Step2-2.GC_correction_ratio.19S6926253.sh
step2_2Dict = dict()
step2_2Dict['name'] = 'Step2-2.GC_correction_ratio'
step2_2Dict['priorStep'] = ['Step2-1.sample']
step2_2Dict['nextStep'] = ['Step2-3']
jobShList = list()
for sampleId in fqpathDict.keys():
	jobDict = dict()
	jobDict['mem'] = '0.5'
	jobDict['sh'] = os.path.join(work_path, 'Shell', 'Step2-2.GC_correction_ratio.{0}.sh'.format(sampleId))
	jobShList.append(jobDict)
step2_2Dict['jobSh'] = jobShList
allJsonList.append(step2_2Dict)

#Step2-3.slope.19S6926253.sh
step2_3Dict = dict()
step2_3Dict['name'] = 'Step2-3.GC_slope'
step2_3Dict['priorStep'] = ['Step2-2.GC_correction_ratio']
step2_3Dict['nextStep'] = ['Step2-4']
jobShList = list()
for sampleId in fqpathDict.keys():
	jobDict = dict()
	jobDict['mem'] = '1'
	jobDict['sh'] = os.path.join(work_path, 'Shell', 'Step2-3.slope.{0}.sh'.format(sampleId))
	jobShList.append(jobDict)
step2_3Dict['jobSh'] = jobShList
allJsonList.append(step2_3Dict)

#Step2-4.19S6926253.sh
step2_4Dict = dict()
step2_4Dict['name'] = 'Step2-4'
step2_4Dict['priorStep'] = ['Step2-3.slope']
step2_4Dict['nextStep'] = ['Step3']
jobShList = list()
for sampleId in fqpathDict.keys():
	jobDict = dict()
	jobDict['mem'] = '1'
	jobDict['sh'] = os.path.join(work_path, 'Shell', 'Step2-4.{0}.sh'.format(sampleId))
	jobShList.append(jobDict)
step2_4Dict['jobSh'] = jobShList
allJsonList.append(step2_4Dict)


# Step3.result.sh
step3Dict = dict()
step3Dict['name'] = 'Step3.result'
step3Dict['priorStep'] = ['Step2']
step3Dict['nextStep'] = ['Step4-1.virus']
jobShList = list()
jobDict = dict()
jobDict['mem'] = '2'
jobDict['sh'] = os.path.join(work_path, 'Shell', 'Step3.result.sh')
jobShList.append(jobDict)
step3Dict['jobSh'] = jobShList
allJsonList.append(step3Dict)

"""
# Step4.annotation.sh
step4Dict = dict()
step4Dict['name'] = 'Step4.Anno'
step4Dict['priorStep'] = ['Step3.result']
step4Dict['nextStep'] = ['Step5.virus']
jobShList = list()
jobDict = dict()
jobDict['mem'] = '2'
jobDict['sh'] = os.path.join(work_path, 'Shell', 'Step4.Anno.sh')
jobShList.append(jobDict)
step4Dict['jobSh'] = jobShList
allJsonList.append(step4Dict)
"""

# virus run
step5Dict = dict()
step5Dict['name'] = 'Step4-1.virus'
step5Dict['priorStep'] = ['Step3.result']
step5Dict['nextStep'] = ['Step4-2.virus.find_positive']
jobShList = list()
for sampleId in fqpathDict.keys():
	jobDict = dict()
	jobDict['mem'] = '6'
	jobDict['sh'] = os.path.join(work_path, 'Shell', 'Step4-1.virus.{0}.sh'.format(sampleId))
	jobShList.append(jobDict)
step5Dict['jobSh'] = jobShList
allJsonList.append(step5Dict)

step6Dict = dict()
step6Dict['name'] = 'Step4-2.virus.find_positive'
step6Dict['priorStep'] = ['Step4-1.virus']
step6Dict['nextStep'] = []
jobShList = list()
jobDict = dict()
jobDict['mem'] = '2'
jobDict['sh'] = os.path.join(work_path, 'Shell', 'Step4-2.virus.find_positive.sh')
jobShList.append(jobDict)
step6Dict['jobSh'] = jobShList
allJsonList.append(step6Dict)


# 没有一体机启动程序，也没有作业调度，只能靠自己启动咯
if auto:
	from auto_qsub import AutoQsub
	autoQsub = AutoQsub()
	autoQsub.start(work_path, allJsonList)

# 格式化json
allJsonList = json.dumps(allJsonList, indent=4)
# print(allJsonList)

allStepsPath = os.path.join(work_path, 'allSteps.json')
with open(allStepsPath, 'w') as f:
	f.write(str(allJsonList))
#logging.info('pscc end')
