# Extract 7,500 random values from the second column of lines starting with 'S_'
#awk '$2 ~ /^S_/{print $2}' /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/3_Linkage/1_chr11/HR.85.maf2.ldwin100k.3500kb.inversion.bim | shuf -n 7500 > 0_selected_values.txt

# Filter HR.ld.gz: remove rows where the 3rd or 6th column matches any value in selected_values.txt
# The filtered content is then re-zipped to HR.ld.filtered.gz
zcat HR.85.maf2.ldwin100k.3500kb.inversion.ld.gz | awk '
BEGIN {
    while (getline < "0_selected_values.txt") {
        values[$1] = 1
    }
}
!(($3 in values) || ($6 in values))' | gzip > HR.85.maf2.ldwin100k.3500kb.inversion.ld.filter.gz

