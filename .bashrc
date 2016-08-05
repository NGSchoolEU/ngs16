
NGSchoolDir=/ngschool

export PATH=$PATH:$NGSchoolDir/bin:$NGSchoolDir/src
export PATH=$PATH:$NGSchoolDir/src/SPAdes-3.9.0-Linux/bin:$NGSchoolDir/src/trinityrnaseq:$NGSchoolDir/src/TransDecoder
export PATH=$PATH:$NGSchoolDir/src/augustus/src::$NGSchoolDir/src/augustus-3.0.2/bin

export PATH=$PATH:$NGSchoolDir/src/snpEff

export PATH=$PATH:$NGSchoolDir/src/Bismark:$NGSchoolDir/src/preseq_v2.0:$NGSchoolDir/src/qualimap_v2.2:$NGSchoolDir/src/trim_galore_zip
export PATH=$PATH::$NGSchoolDir/src/convector

export SSPACEBIN=$NGSchoolDir/src/SSPACE/SSPACE_Standard_v3.0.pl

echo ".bashrc loaded!"
