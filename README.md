# variantFiltering

This script permits variant filtration based on inheritance patterns.




## Usage
```
# git clone https://github.com/mgborges/variantFiltering.git 
cd variantFiltering/test/

vcf="test.vcf"
inheritance="XlinkedRecessive"
phenotype="affected,unaffected,unaffected"
samples="0009709,0015715,0016708"
relation="femaleChild,father,mother"

source ../bin/inheritanceFiltering.sh

filterInheritance $vcf $inheritance $phenotype $samples $relation
```
