cd ..
make clean
make
cd example
sleep 2
ls -l ./ | awk '{print $NF}' | while read frqcase
do
  if [ -d ${frqcase} ]; then
    echo -n "testing ${frqcase} ... "
    start=$(date +%s)
    cd ${frqcase}
    mkdir test
    cp *.dat test/
    cd test
    for folder in DatTemp DatBodyIB DatInfo DatOthe DatBody DatFlow 
    do
      mkdir $folder
    done
    ../../../FSILBM3D &> newrunlog
    cat newrunlog | awk '{if($1=="FIELDSTAT") print $2, $3, $4}' > newfield
    cat ../runlog | awk '{if($1=="FIELDSTAT") print $2, $3, $4}' > oldfield
    result=`diff oldfield newfield`
    end=$(date +%s)
    echo -n "Time used $(($end-$start)) seconds. "
    if [ -z "$result" ]; then
        echo "Passed. "
    else
        echo " Failed: "
        echo $result
    fi
    cd ..
    #rm -rf test
    cd ..
  fi
done
