#!/bin/bash
#---------------------------------------------------------------------------------------------------
# Make CMSSW aware of the qjets external library
#
#
#                                                                   Apr 10, 2015 - V1 Y. Iiyama
#---------------------------------------------------------------------------------------------------
function configureScram {
  EXTERNAL=$1

  XMLPATH=$CMSSW_BASE/config/toolbox/$SCRAM_ARCH/tools/selected/qjets.xml

  # add local qjets external to scram config
  [ -e $XMLPATH ] && mv $XMLPATH ${XMLPATH}-last.$$
  echo \
'
  <tool name="qjets" version="2">
    <info url="http://jets.physics.harvard.edu/Qjets/html/Welcome.html"/>
    <lib name="qjets"/>
    <client>
      <environment name="QJETS_BASE" default="xx-PATH-xx/Qjets"/>
      <environment name="LIBDIR" default="$QJETS_BASE/lib"/>
      <environment name="INCLUDE" default="$QJETS_BASE"/>
    </client>
  </tool>
' | sed "s#xx-PATH-xx#$EXTERNAL#" \
> $XMLPATH

  # commit the scram config changes
  cd $CMSSW_BASE/src
  scram setup qjets
}

# test whether environment is clean and we can start
if [ -z "$CMSSW_BASE" ]
then
  echo ""
  echo " ERROR - cmssw is not setup. Please, setup release directory and checkout MitMonoJet."
  echo ""
  exit 1
fi
if [ -d "$CMSSW_BASE/src/MitMonoJet" ]
then
  echo ""
  echo " INFO - found MitMonoJet location at: $CMSSW_BASE/src/MitMonoJet"
  echo ""
else
  echo ""
  echo " ERROR - MitMonoJet is not in your release. Please check it from GITHUB."
  echo ""
  exit 1
fi

# the default location
EXTERNAL_DEF="/home/cmsprod/cms/external"
if [ -d "$EXTERNAL_DEF" ]
then
  EXTERNAL="$EXTERNAL_DEF"
  echo ""
  echo " INFO - found external location at: $EXTERNAL"
  echo ""
  # use existing location just adjust scram configuration
  configureScram $EXTERNAL
  exit 0
else
  EXTERNAL="/home/$USER/cms/external"
  echo " INFO - default external location ($EXTERNAL_DEF) not found make own external at: $EXTERNAL"
  echo ""
fi

FASTJET=
for ent in $(ls $EXTERNAL_DEF/fastjet-*); do
    if [ -d $ent ]; then
        FASTJET=$ent
        break
    fi
done
if ! [ $FASTJET ]; then
    for ent in $(ls $EXTERNAL/fastjet-*); do
        if [ -d $ent ]; then
            FASTJET=$ent
            break
        fi
    done
fi

if ! [ $FASTJET ]; then
    echo "Qjets is a fastjet plugin. No fastjet package was found in $EXTERNAL_DEF. Aborting."
    exit 1
fi

# Here the real work starts

mkdir -p $EXTERNAL

# download
REMOTEDIR="http://jets.physics.harvard.edu/Qjets"
ZIPNAME="QjetsNew-1"

mkdir -p /tmp/$USER
rm -rf /tmp/$USER/$ZIPNAME 2> /dev/null

echo " INFO - download starting"
echo ""
wget "$REMOTEDIR/${ZIPNAME}.zip" -O /tmp/$USER/${ZIPNAME}.zip

unzip /tmp/$USER/${ZIPNAME}.zip -d /tmp/$USER
rm /tmp/$USER/${ZIPNAME}.zip

rm -rf /tmp/$USER/$ZIPNAME/sample_data 2> /dev/null

mv /tmp/$USER/$ZIPNAME $EXTERNAL/Qjets

# compile
cd $EXTERNAL/Qjets
export PATH=$FASTJET:$PATH
make

# final adjustment to scram configuration
configureScram $EXTERNAL

exit 0
