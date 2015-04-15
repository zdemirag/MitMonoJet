#!/bin/bash

##################################################################
# Move Bambu skim from histogram output directory to /mnt/hadoop #
#                                               Y.I. Apr 13 2015 #
##################################################################

copyRootFile () {
    local src=$1
    local targ=$2
    local current=$3

    local retcode=1
    local file

    for file in $(ls $src/$current)
    do
        if [[ $file =~ .root$ ]]
        then
            echo cp $src/$current/$file $targ/$current/$file
            cp $src/$current/$file $targ/$current/$file
            if [ $? -eq 0 ]
            then
                echo rm $src/$current/$file
                rm $src/$current/$file
            else
                echo "Failed to copy $src/$current/$file"
            fi
            retcode=0
        elif [ -d $src/$current/$file ]
        then
            mkdir -p $targ/$current/$file
            copyRootFile $src $targ $current/$file
            if [ $? -eq 0 ]
            then
                retcode=0
            elif [ $? -eq 1 ]
            then
                # nothing was copied
                rmdir $targ/$current/$file
            fi
        fi
    done

    return $retcode
}


TARGET=/mnt/hscratch/$USER/skim

if [ -e $HOME/cms/skim ]
then
    if [[ $(readlink $HOME/cms/skim) =~ ^/home/$USER ]]
    then
        echo "$HOME/cms/skim cannot be on NFS directory /home/$USER."
        echo "Please move the directory to /mnt/hscratch/$USER."
        exit 1
    fi
else
    mkdir -p $TARGET
    if [ $? -ne 0 ]
    then
        echo "Failed to create directory /mnt/hscratch/$USER/skim."
        exit 1
    fi
    ln -s $TARGET $HOME/cms/skim
fi

if ! [ $MIT_PROD_CFG ] || ! [ $MIT_PROD_BOOK ]
then
    echo "MIT_PROD_CFG or MIT_PROD_BOOK not set"
    exit 1
fi

copyRootFile $HOME/cms/hist/$MIT_PROD_CFG $TARGET/$MIT_PROD_CFG .
