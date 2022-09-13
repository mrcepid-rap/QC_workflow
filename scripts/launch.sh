#!/bin/bash

TAR=$1
PHENO=$2
BIN=$3
OUT=$4
SEX=$5
COVARFILE=$6
QUANTCOVAR=$7
CATCOVAR=$8
TOOL=$9

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

if [[ "$TOOL" = "null" ]]
then
    TOOL=("bolt" "saige" "staar" "glm" "regenie")
else
    TOOL=$(echo $TOOL | tr "," "\n")
fi

DEFAULT_COVARS="--transcript_index file-GFzk5gjJ0zVQXPQX8p4jj2BJ  --base_covariates file-GGPKqYQJJy8Gv2XKJ3v3KK85 --bgen_index file-GFzfyyQJ0zVqXZQBK4592QQg --array_bed_file file-GGJZzg0J80QQgzj2Jk7qX8z4 --array_fam_file file-GGJb0zQJ80QbJ1fgJQpFQpK6 --array_bim_file file-GGJb100J80QVz6627f712Z5G --low_MAC_list file-GGJb138J80QZK57z7fBYjgjv --sparse_grm file-GGJb128J80Qy4fG17p2yKGFq --sparse_grm_sample file-GGJb12jJ80QX0ZJXJ4jxBfGX --inclusion_list file-GGJb118J80QpZGk5Fy4q6bJZ"

if [[ " ${TOOL[*]} " =~ "bolt" ]]
then
    dx run mrcepid-runassociationtesting --brief --yes --name "${OUT}.bolt" --priority normal --destination results/ -imode=burden -ioutput_prefix="${OUT}.bolt" -iinput_args="--association_tarballs $TAR --phenofile $PHENO --sex $SEX --tool bolt --run_marker_tests $APPEND $DEFAULT_COVARS"
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
    dx run mrcepid-runassociationtesting --brief --yes --name "${OUT}.regenie" --priority normal --destination results/ -imode=burden -ioutput_prefix="${OUT}.regenie" -iinput_args="--association_tarballs $TAR --phenofile $PHENO --sex $SEX --tool regenie --run_marker_tests $APPEND $DEFAULT_COVARS"
fi

if [[ " ${TOOL[*]} " =~ "glm" ]]
then
    dx run mrcepid-runassociationtesting --brief --yes --name "${OUT}.glm" --priority normal --destination results/ -imode=burden -ioutput_prefix="${OUT}.glm" -iinput_args="--association_tarballs $TAR --phenofile $PHENO --sex $SEX --tool glm $APPEND $DEFAULT_COVARS"
fi





