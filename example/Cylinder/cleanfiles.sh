for folder in DatTemp DatBodyIB DatBodySpan DatInfo DatOthe DatBody DatFlow 
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

for file in Check.dat log.txt stderr.txt stdout.txt *.*0000
do
  if [[ -e $file ]]; then
    rm $file
    echo Delete $file
  fi
done
