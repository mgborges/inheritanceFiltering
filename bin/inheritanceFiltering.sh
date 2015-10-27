# Split a multi-sample VCF file into separated samples
function splitVCF {
    file=$1
    samples=`less $file | grep "^#CHROM" | cut -f10-`
    samples=($samples)
    numberofSamples=`bc <<< "${#samples[@]} - 1"`
    for i in `seq 0 $numberofSamples`; do
    j=`bc <<< "$i+10"`
    cat $file | cut -f 1-9,$j > /tmp/${samples[$i]}.vcf
    done
}

# Organizes the FORMAT field that contains GT:AD:DP:GQ:PL
function getGenotypes {
    file=$1
    sampleNames=(`cat $file | grep "^#CHROM" | cut -f10-`)
    numberofSamples=${#sampleNames[@]}
    if [ "$numberofSamples" != 1 ]
    then
        echo "Multi-sample file, please split it."
        return 1
    fi
    header="#CHR-POS-REF-ALT"
    echo $header FORMAT $sampleNames GT

# Split the genotype field
    while read -r line
    do
        GT=`echo $line| grep -v ^# | awk '{print $10}' | awk -F: '{print $1}'`
        content=`echo $line| grep -v ^# | awk '{print $1 "-" $2 "-" $4 "-" $5 " " $9 " " $10}'`
        echo $content $GT
    done < $file
}


function getGenotypes_COMPLEX {
    file=$1
    sampleNames=(`cat $file | grep "^#CHROM" | cut -f10-`)
    numberofSamples=${#sampleNames[@]}
    if [ "$numberofSamples" != 1 ]
    then
        echo "Multi-sample file, please split it."
        return 1
    fi

    genotypeAtributes=(`cat $file| grep -v ^# |cut -f9 | sort | uniq | sed 's/:/\n/g' | sort | uniq`)
    numberofgenotypeAtributes=`bc <<< "${#genotypeAtributes[@]} - 1"`

    header=`cat $file| grep ^#CHR | cut -f1,2,4,5 | sed 's/:/\t/g' | sed 's/\t/-/g'`
    echo $header ${genotypeAtributes[*]}

# Split the genotype field
    while read -r line
    do
        if [ `echo $line | grep ^# -c` == "1" ]
            then
            continue
        fi
        sampleColumn="10"
        position=`echo $line| grep -v ^# | cut -d " " -f1,2,4,5 | sed 's/ /_/g'`
        maskField=(`echo $line| grep -v ^# | cut -d " " -f9 | sed 's/:/\t/g'`)
        numberofmaskField=`bc <<< "${#maskField[@]} - 1"`
        genotypeField=(`echo $line| grep -v ^# | cut -d " " -f$sampleColumn | sed 's/:/\t/g'`)

        fieldInfo=". . . . . . ."
        fieldInfo=($fieldInfo)

        for a in `seq 0 $numberofgenotypeAtributes`
            do
                 for m in `seq 0 $numberofmaskField`
                     do
                      if [ "${genotypeAtributes[$a]}" = "${maskField[$m]}" ] && [ "${genotypeField[0]}" != "./." ]
                      then
                           fieldInfo[$a]=${genotypeField[$m]}
                      fi
                     done
            done

        echo $position ${fieldInfo[*]}
    done < $file
}




# Which is the FORMAT colum that contains a certain value?
function whichColum
{
     filename=$1
     field=$2
     header=`cat $filename | grep ^#`
     header=($header)
     fieldNumber=""

     for i in `seq 1 ${#header[@]}`
     do
          let i=$i-1
          if [ ${header[$i]} = $field ]
          then
               let fieldNumber=$i+1
          fi
     done
     echo $fieldNumber
}

# Homozygous variants 0/0
function isHomo
{
     grep "0/0" $1 | cut -f 1 -d " "
}

# Heterozygous variants 0/1 1/2 ...
function isHetero
{
     columnNumber=`whichColum $1 GT`     
     while read -r line
     do
     first=`echo $line | cut -d " " -f $columnNumber | awk -F '[/|]' '{print $1}'`
     second=`echo $line | cut -d " " -f $columnNumber | awk -F '[/|]' '{print $2}'`
     if [ "$first" != "$second" ]
     then
          echo $line | grep "/\||" -m 1 | cut -f 1 -d " "
     fi
     done < $1 
}


# Homozygous alternative variants 1/1 2/2 ...
function isHomoAlter
{
     columnNumber=`whichColum $1 GT`
     while read -r line
     do
     first=`echo $line | cut -d " " -f $columnNumber | awk -F '[/|]' '{print $1}'`
     second=`echo $line | cut -d " " -f $columnNumber | awk -F '[/|]' '{print $2}'`
     if [ "$first" = "$second" ]
     then
          echo $line | grep -v "0/0" | grep "/\||" -m 1 | cut -f 1 -d " "
     fi
     done < $1 
}

# Given a *.phenotype file, it splits into the possible genotypes
function splitPhenotypes
{
     isHomo $1 > $1.homo
     isHetero $1 > $1.hetero
     isHomoAlter $1 > $1.homoAlter
}

# Union of two files
function union
{
     file1=$1
     file2=$2
     cat $file1 $file2 | sort | uniq
}

# Intersection of two files
function intersection
{
     file1=$1
     file2=$2
     fgrep -wf $file1 $file2
}

# Difference of files
function difference
# lines in FILE1 that are not in FILE2
{
     file1=$1
     file2=$2
     fgrep -wvf $file2 $file1
}

# AUTOSSOMAL DOMINANT INHERITANCE
function autossomalDominant
{
phenotype=$1
phenotype=`echo $phenotype | sed 's/,/ /g'`
phenotype=($phenotype)
phenotypeComparations=`bc <<< "${#phenotype[@]} - 1"`

samples=$2
samples=`echo $samples | sed 's/,/ /g'`
samples=($samples)

sample1=${samples[0]}

file1=`ls /tmp/*${samples[0]}*.hetero`
anteriorPhenotype=${phenotype[0]}

for i in `seq 1 $phenotypeComparations`
do
     file2=`ls /tmp/*${samples[$i]}*.hetero`
     output=/tmp/$i.tmp
     if [ $anteriorPhenotype == "affected" ] && [ ${phenotype[$i]} == "affected" ]
     then
          intersection $file1 $file2 > $output
          anteriorPhenotype="affected"
     elif [ $anteriorPhenotype == "affected" ] && [ ${phenotype[$i]} == "unaffected" ]
     then
          difference $file1 $file2 > $output
          anteriorPhenotype="affected"
     elif [ $anteriorPhenotype == "unaffected" ] && [ ${phenotype[$i]} == "unaffected" ]
     then
          union $file1 $file2 > $output
          anteriorPhenotype="unaffected"
     elif [ $anteriorPhenotype == "unaffected" ] && [ ${phenotype[$i]} == "affected" ]
     then
          difference $file2 $file1 > $output
          anteriorPhenotype="affected"
     fi
     file1=$output
done
cat $file1
}

# Autossomal recessive inheritance
function autossomalRecessive
{
phenotype=$1
phenotype=`echo $phenotype | sed 's/,/ /g'`
phenotype=($phenotype)
phenotypeComparations=`bc <<< "${#phenotype[@]} - 1"`

samples=$2
samples=`echo $samples | sed 's/,/ /g'`
samples=($samples)

sample1=${samples[0]}

file1=`ls /tmp/*${samples[0]}*.homoAlter`
anteriorPhenotype=${phenotype[0]}

for i in `seq 1 $phenotypeComparations`
do
     file2=`ls /tmp/*${samples[$i]}*.homoAlter`
     output=/tmp/$i.tmp
     if [ $anteriorPhenotype == "affected" ] && [ ${phenotype[$i]} == "affected" ]
     then
          intersection $file1 $file2 > $output
          anteriorPhenotype="affected"
     elif [ $anteriorPhenotype == "affected" ] && [ ${phenotype[$i]} == "unaffected" ]
     then
          difference $file1 $file2 > $output
          anteriorPhenotype="affected"
     elif [ $anteriorPhenotype == "unaffected" ] && [ ${phenotype[$i]} == "unaffected" ]
     then
          union $file1 $file2 > $output
          anteriorPhenotype="unaffected"
     elif [ $anteriorPhenotype == "unaffected" ] && [ ${phenotype[$i]} == "affected" ]
     then
          difference $file2 $file1 > $output
          anteriorPhenotype="affected"
     fi
     file1=$output
done
cat $file1
}

# Mitochondrial inheritance
function mito
{
phenotype=$1
phenotype=`echo $phenotype | sed 's/,/ /g'`
phenotype=($phenotype)
phenotypeComparations=`bc <<< "${#phenotype[@]} - 1"`

samples=$2
samples=`echo $samples | sed 's/,/ /g'`
samples=($samples)

sample1=${samples[0]}

file1=`ls /tmp/*${samples[0]}*.homoAlter`
anteriorPhenotype=${phenotype[0]}

for i in `seq 1 $phenotypeComparations`
do
     file2=`ls /tmp/*${samples[$i]}*.homoAlter`
     output=/tmp/$i.tmp
     if [ $anteriorPhenotype == "affected" ] && [ ${phenotype[$i]} == "affected" ]
     then
          intersection $file1 $file2 > $output
          anteriorPhenotype="affected"
     elif [ $anteriorPhenotype == "affected" ] && [ ${phenotype[$i]} == "unaffected" ]
     then
          difference $file1 $file2 > $output
          anteriorPhenotype="affected"
     elif [ $anteriorPhenotype == "unaffected" ] && [ ${phenotype[$i]} == "unaffected" ]
     then
          union $file1 $file2 > $output
          anteriorPhenotype="unaffected"
     elif [ $anteriorPhenotype == "unaffected" ] && [ ${phenotype[$i]} == "affected" ]
     then
          difference $file2 $file1 > $output
          anteriorPhenotype="affected"
     fi
     file1=$output
done
cat $file1 | grep ^MT
}


# XLINKED DOMINANT INHERITANCE
# FATHER AFFECTED
function XlinkedDominantFatherAffected
{
phenotype=$1
phenotype=`echo $phenotype | sed 's/,/ /g'`
phenotype=($phenotype)
phenotypeComparations=`bc <<< "${#phenotype[@]} - 1"`

samples=$2
samples=`echo $samples | sed 's/,/ /g'`
samples=($samples)

relation=$3
relation=`echo $relation | sed 's/,/ /g'`
relation=($relation)
relationNumber=`bc <<< "${#relation[@]} - 1"`

files=""

for i in `seq 0 $relationNumber`
do
     r=${relation[$i]}
     if [ "$r" == "femaleChild" ] && [ ${phenotype[$i]} == "affected" ]
     then
          file=`ls /tmp/*${samples[$i]}*.hetero`
          files=`echo $files $file`
     elif [ "$r" == "father" ] && [ ${phenotype[$i]} == "affected" ]
     then
          file=`ls /tmp/*${samples[$i]}*.homoAlter`
          files=`echo $files $file`
     elif [ "$r" == "mother" ] && [ ${phenotype[$i]} == "unaffected" ]
     then
          file=`ls /tmp/*${samples[$i]}*.homo`
          files=`echo $files $file`
     else
          echo ERROR $r should not be ${phenotype[$i]}; return 1
     fi
let i=$i+1
done

files=($files)

file1=${files[0]}

for i in `seq 1 $phenotypeComparations`
do
     file2=${files[$i]}
     output=/tmp/$i.tmp
     intersection $file1 $file2 > $output
     file1=$output
done
cat $file1 | grep ^X
}

# Returns if the affected child is male or female
function whoIsTheAffectedChild
{
relation=$1
phenotype=$2
relationNumber=`bc <<< "${#relation[@]} - 1"`
for i in `seq 0 $relationNumber`
do
     r=${relation[$i]}
     if [ "$r" == "femaleChild" ] && [ ${phenotype[$i]} == "affected" ]
     then 
         echo 0; return 0
     elif [ "$r" == "maleChild" ] && [ ${phenotype[$i]} == "affected" ]
     then
         echo 1; return 0
     fi
done
}

# XLINKED DOMINANT INHERITANCE
# MOTHER AFFECTED
function XlinkedDominantMotherAffected
{
phenotype=$1
phenotype=`echo $phenotype | sed 's/,/ /g'`
phenotype=($phenotype)
phenotypeComparations=`bc <<< "${#phenotype[@]} - 1"`

samples=$2
samples=`echo $samples | sed 's/,/ /g'`
samples=($samples)

relation=$3
relation=`echo $relation | sed 's/,/ /g'`
relation=($relation)
relationNumber=`bc <<< "${#relation[@]} - 1"`

files=""

familialProfile=`whoIsTheAffectedChild $relation $phenotype`

if [ "$familialProfile" == 0 ] # femaleChild AFFECTED
then
     for i in `seq 0 $relationNumber`
     do
          r=${relation[$i]}
          if [ "$r" == "femaleChild" ] && [ ${phenotype[$i]} == "affected" ]
          then
               file=`ls /tmp/*${samples[$i]}*.hetero`
               files=`echo $files $file`
          elif [ "$r" == "father" ] && [ ${phenotype[$i]} == "unaffected" ]
          then
               file=`ls /tmp/*${samples[$i]}*.homo`
                    files=`echo $files $file`
          elif [ "$r" == "mother" ] && [ ${phenotype[$i]} == "affected" ]
          then
               file=`ls /tmp/*${samples[$i]}*.hetero`
               files=`echo $files $file`
          else
               echo ERROR $r should not be ${phenotype[$i]}; return 1
          fi
     let i=$i+1
     done
elif [ "$familialProfile" == 1 ] # maleChild AFFECTED
then
     for i in `seq 0 $relationNumber`
     do
          r=${relation[$i]}
          if [ "$r" == "maleChild" ] && [ ${phenotype[$i]} == "affected" ]
          then
               file=`ls /tmp/*${samples[$i]}*.homoAlter`
               files=`echo $files $file`
          elif [ "$r" == "father" ] && [ ${phenotype[$i]} == "unaffected" ]
          then
               file=`ls /tmp/*${samples[$i]}*.homo`
                    files=`echo $files $file`
          elif [ "$r" == "mother" ] && [ ${phenotype[$i]} == "affected" ]
          then
               file=`ls /tmp/*${samples[$i]}*.hetero`
               files=`echo $files $file`
          else
               echo ERROR $r should not be ${phenotype[$i]}; return 1
          fi
     let i=$i+1
     done
else
     return 1
fi

files=($files)

file1=${files[0]}

for i in `seq 1 $phenotypeComparations`
do
     file2=${files[$i]}
     output=/tmp/$i.tmp
     intersection $file1 $file2 > $output
     file1=$output
done
cat $file1 | grep ^X
}

# Returns if the father is affected or not
function IsTheFatherAffected
{
relation=$1
phenotype=$2
relationNumber=`bc <<< "${#relation[@]} - 1"`
for i in `seq 0 $relationNumber`
do
     r=${relation[$i]}
     if [ "$r" == "father" ] && [ ${phenotype[$i]} == "affected" ]
     then 
         echo 1; return 0
     elif [ "$r" == "father" ] && [ ${phenotype[$i]} != "affected" ]
     then
         echo 0; return 0
     fi
done
}

# X linked recessive inheritance
function XlinkedRecessive
{
phenotype=$1
phenotype=`echo $phenotype | sed 's/,/ /g'`
phenotype=($phenotype)
phenotypeComparations=`bc <<< "${#phenotype[@]} - 1"`

samples=$2
samples=`echo $samples | sed 's/,/ /g'`
samples=($samples)

relation=$3
relation=`echo $relation | sed 's/,/ /g'`
relation=($relation)
relationNumber=`bc <<< "${#relation[@]} - 1"`

files=""

familialProfile=`IsTheFatherAffected $relation $phenotype`

if [ "$familialProfile" == "0" ]
then
     sample1=${samples[0]}

     file1=`ls /tmp/*${samples[0]}*.homoAlter`
     anteriorPhenotype=${phenotype[0]}

     for i in `seq 1 $phenotypeComparations`
     do
          file2=`ls /tmp/*${samples[$i]}*.homoAlter`
          output=/tmp/$i.tmp
          if [ $anteriorPhenotype == "affected" ] && [ ${phenotype[$i]} == "affected" ]
          then
               intersection $file1 $file2 > $output
               anteriorPhenotype="affected"
          elif [ $anteriorPhenotype == "affected" ] && [ ${phenotype[$i]} == "unaffected" ]
          then
               difference $file1 $file2 > $output
               anteriorPhenotype="affected"
          elif [ $anteriorPhenotype == "unaffected" ] && [ ${phenotype[$i]} == "unaffected" ]
          then
               union $file1 $file2 > $output
               anteriorPhenotype="unaffected"
          elif [ $anteriorPhenotype == "unaffected" ] && [ ${phenotype[$i]} == "affected" ]
          then
               difference $file2 $file1 > $output
               anteriorPhenotype="affected"
          fi
          file1=$output
     done
elif [ "$familialProfile" == "1" ]
then
     for i in `seq 0 $relationNumber`
     do
          r=${relation[$i]}
          if [ "$r" == "femaleChild" ]
          then
               file=`ls /tmp/*${samples[$i]}*.hetero`
               files=`echo $files $file`
          elif [ "$r" == "father" ] && [ ${phenotype[$i]} == "affected" ]
          then
               file=`ls /tmp/*${samples[$i]}*.homoAlter`
               files=`echo $files $file`
          else
               file=`ls /tmp/*${samples[$i]}*.homo`
               files=`echo $files $file`
          fi
     let i=$i+1
     done

     files=($files)

     file1=${files[0]}

     for i in `seq 1 $phenotypeComparations`
     do
          file2=${files[$i]}
          output=/tmp/$i.tmp
          intersection $file1 $file2 > $output
          file1=$output
     done
fi
cat $file1 | grep ^X
}

# Given the variants, returns the VCF file
function getVCF
{
vcf=$1
grep ^# $vcf
while read -r line
do
     chr=`echo $line | awk -F - '{print $1}'`
     pos=`echo $line | awk -F - '{print $2}'`
     ref=`echo $line | awk -F - '{print $3}'`
     alt=`echo $line | awk -F - '{print $4}'`
     grep -w ^$chr $vcf | grep -w $pos | grep -w $ref | grep -w $alt
done < $2
}

##################################
### MAIN #########################
##################################
function filterInheritance
{
vcf=`echo $@ | awk -F " " '{print $1}'`
inheritance=`echo $@ | awk -F " " '{print $2}'`
phenotype=`echo $@ | awk -F " " '{print $3}'`
SAMPLES="`echo $@ | awk -F " " '{print $4}'`"
relation=`echo $@ | awk -F " " '{print $5}'`

isThereHelp=`echo $@ | grep '\-help' -c`

if [ "$isThereHelp" == "1" ]
then
echo "filterInheritance filters inheritance profiles in samples of a VCF file."
echo 
echo "USAGE: filterInheritance <VCF> <inheritance> <phenotype> <samples> <relation>"
echo 
echo "VCF - is a multisample VCF file"
echo "inheritance - the type of filter to apply. (autossomalDominant, autossomalRecessive, mito, XlinkedDominantFatherAffected, XlinkedDominantMotherAffected, XlinkedRecessive)"
echo "phenotype - comma separated list of phenotypes (affected or unaffected)"
echo "relation - comma separated list of realations between samples (femaleChild, maleChild, father, mother)"
echo "E.G.:"
echo "filterInheritance FILE.vcf XlinkedRecessive affected,unaffected,unaffected 0009709,0015715,0016708 femaleChild,father,mother"
fi

splitVCF $vcf

vcfs=""
splitSamples=`echo $SAMPLES | sed 's/,/ /g'`

for i in $splitSamples
do
v=`ls /tmp/*$i*.vcf`
vcfs=`echo $vcfs $v`
done

for v in $vcfs
do
getGenotypes $v > $v.phenotype
splitPhenotypes $v.phenotype
done
$inheritance $phenotype $SAMPLES $relation > /tmp/$vcf.pos.tmp
getVCF $vcf /tmp/$vcf.pos.tmp
}
#### end ###
