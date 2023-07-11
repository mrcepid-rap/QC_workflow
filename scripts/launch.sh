#!/bin/bash

TAR=$1
PHENO=$2
PHENONAME=$3
BIN=$4
OUT=$5
SEX=$6
COVARFILE=$7
QUANTCOVAR=$8
CATCOVAR=$9
TOOL=${10}
INCLUSION=${11}
EXCLUSION=${12}


if [[ "$BIN" == "true" ]]
then
    APPEND="--is_binary"
elif [[ "$BIN" == "false" ]]
then
     APPEND=""
else
    echo "is_binary input not true/false."
    exit 1
fi

if [[ "$COVARFILE" = "null" ]]
then
    APPEND="$APPEND"
else
    if [[ "$QUANTCOVAR" != "null" ]] && [[ "$CATCOVAR" != "null" ]]
    then
	APPEND="$APPEND --covarfile=$COVARFILE --quantitative_covariates $QUANTCOVAR --categorical_covariates $CATCOVAR"
    elif [[ "$QUANTCOVAR" != "null" ]] && [[ "$CATCOVAR" = "null" ]]
    then
	APPEND="--covarfile $COVARFILE --quantitative_covariates $QUANTCOVAR"
    elif [[ "$QUANTCOVAR" = "null" ]] && [[ "$CATCOVAR" != "null" ]]
    then
	APPEND="$APPEND --covarfile $COVARFILE --categorical_covariates $CATCOVAR"
    else
	echo "Supplied a covariates file but didn't specify covariate names."
	exit 1
	APPEND=""
    fi
fi

if [[ "$PHENONAME" != "null" ]]
then
   APPEND="$APPEND --phenoname $PHENONAME"
fi

if [[ "$INCLUSION" != "null" ]]
then
   APPEND="$APPEND --inclusion_list $INCLUSION"
fi

if [[ "$EXCLUSION" != "null" ]]
then
   APPEND="$APPEND --exclusion_list $EXCLUSION"
fi

if [[ "$TOOL" = "null" ]]
then
    TOOL=("bolt" "saige" "staar" "glm" "regenie")
else
    TOOL=$(echo $TOOL | tr "," "\n")
fi

DEFAULT_COVARS="--transcript_index file-GFzk5gjJ0zVQXPQX8p4jj2BJ --base_covariates file-GGZkYk8JJy8GFjjF6kYG01g8 --bgen_index file-GFzfyyQJ0zVqXZQBK4592QQg --array_bed_file file-GGbZKfjJP7J5223J4v6k49j0 --array_fam_file file-GGbZPy0JP7JJBjkG4vbPqKk9 --array_bim_file file-GGbZPyQJP7J82gv54x6F3Kjq --low_MAC_list file-GGbZQ08JP7JJPBbP4vfF9X25 --sparse_grm file-GGbZPz8JP7J4vVXp4vj5JGyx --sparse_grm_sample file-GGbZPzjJP7JF98v34vxg8QkV"

if [[ " ${TOOL[*]} " =~ "bolt" ]]
then
    dx run mrcepid-runassociationtesting --brief --yes --name "${OUT}.bolt" --priority normal --destination results/ -imode=burden -ioutput_prefix="${OUT}.bolt" -iinput_args="--association_tarballs $TAR --phenofile $PHENO --sex $SEX --tool bolt --run_marker_tests --bolt_non_infinite $APPEND $DEFAULT_COVARS"
fi

if [[ " ${TOOL[*]} " =~ "saige" ]]
then
    dx run mrcepid-runassociationtesting --brief --yes --name "${OUT}.saige" --priority normal --destination results/ -imode=burden -ioutput_prefix="${OUT}.saige" -iinput_args="--association_tarballs $TAR --phenofile $PHENO --sex $SEX --tool saige $APPEND $DEFAULT_COVARS"
fi

if [[ " ${TOOL[*]} " =~ "staar" ]]
then
    dx run mrcepid-runassociationtesting --brief --yes --name "${OUT}.staar" --priority normal --destination results/ -imode=burden -ioutput_prefix="${OUT}.staar" -iinput_args="--association_tarballs $TAR --phenofile $PHENO --sex $SEX --tool staar $APPEND $DEFAULT_COVARS"
fi

if [[ " ${TOOL[*]} " =~ "regenie" ]]
then
    dx run mrcepid-runassociationtesting --brief --yes --name "${OUT}.regenie" --priority high --destination results/ -imode=burden -ioutput_prefix="${OUT}.regenie" -iinput_args="--association_tarballs $TAR --phenofile $PHENO --sex $SEX --tool regenie --run_marker_tests $APPEND $DEFAULT_COVARS"
fi

if [[ " ${TOOL[*]} " =~ "glm" ]]
then
    dx run mrcepid-runassociationtesting --brief --yes --name "${OUT}.glm" --priority normal --destination results/ -imode=burden -ioutput_prefix="${OUT}.glm" -iinput_args="--association_tarballs $TAR --phenofile $PHENO --sex $SEX --tool glm $APPEND $DEFAULT_COVARS"
fi





