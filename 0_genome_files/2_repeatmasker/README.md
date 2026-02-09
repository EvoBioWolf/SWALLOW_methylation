## Repeat Filtering

I wanted to remove repeat regions, as these are harder to map with short read data and therefore, I wanted to avoid mis-mapping in those regions.

I identified repeats using 0_repeat_masker.sh

I then use 1_gff2bed.sh to processes the **GFF** (General Feature Format) file, which contains the repeat annotations to create a simplified **BED** file.