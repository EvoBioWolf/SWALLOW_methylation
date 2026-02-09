gunzip -c /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/3_Linkage/1_Rall/HR.Rall.maf2.ldwin100.70kb.ld.gz | awk '{
    if ($3 ~ /^M_/ && $6 ~ /^M_/) {
        print > "HR_Rall_meth.maf2.ldwin100"
    } else if ($3 ~ /^S_/ && $6 ~ /^S_/) {
        print > "HR_Rall_snp.maf2.ldwin100"
    } else {
        print > "HR_Rall_snp.meth.maf2.ldwin100"
    }
}' 

gzip HR_Rall_meth.maf2.ldwin100
gzip HR_Rall_snp.maf2.ldwin100
gzip HR_Rall_snp.meth.maf2.ldwin100
