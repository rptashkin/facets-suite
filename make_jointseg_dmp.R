#!/usr/bin/env Rscript

pkgs = c('data.table', 'argparse')
lapply(pkgs, require, character.only = T)

parser=ArgumentParser(python_cmd="/dmp/resources/prod/tools/system/python/production/bin/python")  
parser$add_argument('-i', '--id', type='character', help='sampleID')
parser$add_argument('-w', '--wd', type='character', help='WorkingDir')
args=parser$parse_args()
  
id <- args$id
wd <- args$wd

outfile <- paste(id,"hisens.Rdata.jointseg", sep="_")
outfile <- paste(wd,outfile, sep="/")

infile = paste(id,"hisens.Rdata", sep="_")
infile = paste(wd,infile, sep="/")

load(infile)
write.table(out$jointseg,file=outfile,sep="\t",row.names=FALSE)
