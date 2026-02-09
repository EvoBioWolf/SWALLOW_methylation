gunzip -c /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/3_Linkage/1_chr11/HR.85.maf2.ldwin100.70kb.ld.gz | awk '{
    if ($3 ~ /^M_/ && $6 ~ /^M_/) {
        print > "HR_85_meth.maf2.ldwin100"
    } else if ($3 ~ /^S_/ && $6 ~ /^S_/) {
        print > "HR_85_snp.maf2.ldwin100.t"
    } else {
        print > "HR_85_snp.meth.maf2.ldwin100"
    }
}' 

gzip HR_85_meth.maf2.ldwin100
rm HR_85_snp.maf2.ldwin100.t
gzip HR_85_snp.meth.maf2.ldwin100
