if [ -z $CMSSW_BASE ]
then
  echo ""
  echo " Setting up MitHgg failed! (\$CMSSW_BASE = is empty)."
  echo ""
else
  export MIT_MONOJET_DIR="$CMSSW_BASE/src/MitMonoJet"
  export PATH="$MIT_MONOJET_DIR/bin:${PATH}"
  export PYTHONPATH="$MIT_MONOJET_DIR/python:${PYTHONPATH}"
fi

