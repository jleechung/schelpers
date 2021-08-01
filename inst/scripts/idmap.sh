
################################################################
##### Help #####################################################
################################################################

help()
{
	echo "Fetch mapping between Ensembl gene ID, Ensembl transcript ID and common gene name."
	echo
	echo "Usage: idmap [-i|o]"
	echo "Options:"
	echo "i		Input gtf file"
	echo "o		Output map file (default: map.txt)"
}

tlog()
{
	date +"%T"
}


################################################################
##### Parse arguments ##########################################
################################################################

output="map.txt"

while getopts ":hi:o:" option; do
	case $option in
		h)
			help
			exit;;
		i)
			input=$OPTARG;;
		o)
			output=$OPTARG;;
		\?)
			echo "Invalid option. -h for usage"
			exit;;
	esac
done

if [[ $# -eq 0 ]]
then
	help
	exit 1
fi

################################################################
##### Main #####################################################
################################################################

## Check the input file
if [[ ! -e "$input" ]]
then
	echo "Error: ${input} does not exist..."
	exit 1
fi

if [[ $input != *gtf.gz ]] && [[ $input != *.gtf ]]
then
	echo "Error: ${input} does not seem to be a gtf file..."
	exit 1
fi

## Parse input

echo "[$(tlog)] Parsing ${input}..."

if [[ $input == *gz ]]
then
	zcat $input | grep -v '^#' | awk -F'\t' '{if ( $3 == "transcript") print}' | \
				cut -f 9 | cut -d' ' -f 2,4,8 | tr -d '"|;' > $output
else
	cat $input | grep -v '^#' | awk -F'\t' '{if ( $3 == "transcript") print}' | \
		cut -f 9 | cut -d' ' -f 2,4,8 | tr -d '"|;' > $output
fi

echo "[$(tlog)] Writing results to ${output}..."

exit 0


