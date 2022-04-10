for folder in DatTemp DatBodyIB DatInfo DatOthe DatBody DatFlow 
do
  echo $folder
  if [[ -d $folder ]]; then
    rm $folder/*
  else
    mkdir $folder
  fi
done
