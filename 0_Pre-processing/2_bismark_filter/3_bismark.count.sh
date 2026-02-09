#Run within /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2021_SwallowRRBS/0_Pre-processing/2_bismark_filter/1_output_files
#Then output can be plotted in R

for i in *COUNTS.txt; do tr "\n" "\t" < $i > $i".sorted"; done
rm *COUNTS.txt
awk '{print FILENAME" \"" $0"\""; nextfile}' *COUNTS* >> a.COUNTS.merged.tmp
sed -i 's/"//g' a.COUNTS.merged.tmp
awk '{print $1, $3, $5, $7, $9, $11, $13, $15}' a.COUNTS.merged.tmp >> a.COUNTS.clean.txt
sed -i 's/.COUNTS.txt.sorted//g' a.COUNTS.clean.txt
sed -i 's/__.................................................//'  a.COUNTS.clean.txt

rm *tmp*