#!/bin/bash
set -e
for files in `ls -d IMG*`;do
        cd $files
        echo $files
		subfile=`ls -d IMG_Data/*assembled.fna`
		subpath=`ls -d IMG_Data/*/`
		echo $subpath
		tag=${subfile##*/}
		tag=${tag%.assembled.fna}
		echo $subfile
		# grep ">" $subfile | sed "s/>//g" | awk '{print $1, $2}' > $subpath/${tag}.gene_oid_2_seq_id.txt
		# grep ">" $subfile | sed "s/>//g" | awk -v tag2=$tag '{print tag2}' > $subpath/${tag}.contig_id.txt
		cat $subpath/${tag}.contig_id.txt >> ../combined_contig_labels.txt
        cd -
done
