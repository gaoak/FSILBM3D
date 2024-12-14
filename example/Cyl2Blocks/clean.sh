
    for folder in DatTemp DatBodyIB DatBodySpan DatInfo DatOthe DatBody DatFlow 
    do
      if [ -d $folder ]; then
        cd $folder
        rm -rf *
        cd ..
      else
        mkdir $folder
      fi
    done
