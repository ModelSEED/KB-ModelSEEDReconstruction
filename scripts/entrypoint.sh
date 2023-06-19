#!/bin/bash

. /kb/deployment/user-env.sh

python ./scripts/prepare_deploy_cfg.py ./deploy.cfg ./work/config.properties

if [ -f ./work/token ] ; then
  export KB_AUTH_TOKEN=$(<./work/token)
fi

if [ $# -eq 0 ] ; then
  sh ./scripts/start_server.sh
elif [ "${1}" = "test" ] ; then
  echo "Run Tests"
  make test
elif [ "${1}" = "async" ] ; then
  sh ./scripts/run_async.sh
elif [ "${1}" = "init" ] ; then
  echo "Initialize module"
  cd /data
  curl -L https://bioseed.mcs.anl.gov/~fliu/modelseedpy/knn_ACNP_RAST_full_01_17_2023_features.json > knn_ACNP_RAST_full_01_17_2023_features.json
  curl -L https://bioseed.mcs.anl.gov/~fliu/modelseedpy/knn_ACNP_RAST_full_01_17_2023.pickle > knn_ACNP_RAST_full_01_17_2023.pickle 
 
  if [ -f knn_ACNP_RAST_full_01_17_2023_features.json ] && [ -f knn_ACNP_RAST_full_01_17_2023.pickle ]; then
  	touch __READY__
  else
    echo "Init failed"
  fi
elif [ "${1}" = "bash" ] ; then
  bash
elif [ "${1}" = "report" ] ; then
  export KB_SDK_COMPILE_REPORT_FILE=./work/compile_report.json
  make compile
else
  echo Unknown
fi
