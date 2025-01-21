for folder in DatBodySpan DatInfo DatBody DatFlow 
do
  if [[ -d $folder ]]; then
    rm -r $folder
    mkdir $folder
    echo Clear $folder
  else
    mkdir $folder
    echo Create $folder
  fi
done

for file in check.dat continue.dat log.* *.txt
do
  if [[ -e $file ]]; then
    rm $file
    echo Delete $file
  fi
done
