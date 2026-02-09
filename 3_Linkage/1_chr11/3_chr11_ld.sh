gunzip -c /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/3_Linkage/1_Rall/HR_Rall_snp.meth.maf2.ldwin100.gz | grep 'chr11' > 3_chr11_maf2.ldwin100.70kb.tmp.txt

awk '$2 >= 1.8e7 && $2 <= 2.0e7' 3_chr11_maf2.ldwin100.70kb.tmp.txt > 3_chr11_maf2.ldwin100.70kb.txt

rm *tmp*