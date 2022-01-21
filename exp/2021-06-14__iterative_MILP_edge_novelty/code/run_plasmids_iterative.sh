for id in {11..60}
do
	echo $id
	time python links_to_fasta.py $id
done	
