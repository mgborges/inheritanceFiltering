# variantFiltering

This script permits variant filtration based on inheritance patterns.

## Usage

`filterInheritance <VCF> <inheritance> <phenotype> <samples> <relation>`

`<VCF>` - is a multisample VCF file

`<inheritance>` - the type of filter to apply:

* `autossomalDominant`;
* `autossomalRecessive`;
* `mito`;
* `XlinkedDominantFatherAffected`;
* `XlinkedDominantMotherAffected`;
* `XlinkedRecessive`.

`<phenotype>` - comma separated list of phenotypes (`affected` or `unaffected`)

`<relation>` - comma separated list of realations between samples (`femaleChild`, `maleChild`, `father`, `mother`)"

## Example
```
# git clone https://github.com/mgborges/variantFiltering.git 
cd variantFiltering/test/

vcf="test.vcf"
inheritance="XlinkedRecessive"
phenotype="affected,unaffected,unaffected"
samples="0009709,0015715,0016708"
relation="femaleChild,father,mother"

source ../bin/inheritanceFiltering.sh

# takes ~10 minutes
filterInheritance $vcf $inheritance $phenotype $samples $relation
```
