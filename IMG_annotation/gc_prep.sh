#!/bin/bash
set -e
for files in `ls -d IMG*`;do
        cd $files
        echo $files
		subfile=`ls -d IMG_Data/*/*genes.faa`
		subpath=`ls -d IMG_Data/*/`
		tag=${subfile##*/}
		tag=${tag%.genes.faa}
		grep ">" $subfile | sed "s/>//g" | awk '{print $1, $2}' > $subpath/${tag}.gene_oid_2_seq_id.txt
        cd -
done
