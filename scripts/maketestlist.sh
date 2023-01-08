#!bin/bash

nfiles=${1:-10}

rm -f runs_local_test.list

#get_file_list.pl -cond trgsetupname=production_pAu200_2015,filename~st_physics,filetype=daq_reco_picoDst,storage=local,production=P16id -keys 'path,filename' -delim / -limit 0 -o runs_path_all.list
#get_file_list.pl -keys path,filename -cond trgsetupname=production_pAu200_2015,library=SL20d,production=P18ih,filetype=daq_reco_picoDst,storage!=HPSS -limit 105 -delim "/" | sed 's|^/home/starlib|root://xrdstar.rcf.bnl.gov:1095//home/starlib|g' > tmp.list


get_file_list.pl -keys 'path,filename' -cond trgsetupname=production_pAu200_2015,filename~st_physics,production=P16id,filetype=daq_reco_picoDst,storage!=HPSS -limit 105 -delim "/" | sed 's|^/home/starlib|root://xrdstar.rcf.bnl.gov:1095//home/starlib|g' > tmp.list



tail -n $nfiles tmp.list > runs_local_test.list
rm -f tmp.list 

