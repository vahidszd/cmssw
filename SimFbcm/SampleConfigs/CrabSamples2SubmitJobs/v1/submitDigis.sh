eval `scramv1 runtime -sh`;
source /cvmfs/cms.cern.ch/common/crab-setup.sh;
for pu in 0p5
do 
    crab submit -c crabDIGI_cfg.py General.requestName=FBCMDIGIPU$pu

done
