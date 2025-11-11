#interactive node
salloc -A gts-ktk3-coda20 -qembers -t 3:00:00 --ntasks-per-node=4

sed 's/\t/,/g' file.tsv > file.csv



for file in *.tsv; do
  sed 's/\t/,/g' "$file" > "${file%.tsv}.csv"
done
