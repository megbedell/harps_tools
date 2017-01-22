for i in *.tar; do mv "$i" "`echo $i | sed -e 's,:,-,g'`"; done
ls *.tar | xargs -n1 tar -xf
find data/* -type f -print | xargs -I{} mv {} .
rm -r data
rm ADP*