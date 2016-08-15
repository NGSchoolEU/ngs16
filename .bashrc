
NGSchoolDir=/ngschool

export PATH=$NGSchoolDir/src/FastQC:$PATH:$NGSchoolDir/src/bin:$NGSchoolDir/src
export PATH=$PATH:$NGSchoolDir/src/SPAdes-3.9.0-Linux/bin:$NGSchoolDir/src/quast
export PATH=$PATH:$NGSchoolDir/src/trinityrnaseq:$NGSchoolDir/src/TransDecoder
export PATH=$PATH:$NGSchoolDir/src/STAR/bin/Linux_x86_64:$NGSchoolDir/src/salmon/build/src:$NGSchoolDir/src/cufflinks-2.2.1.Linux_x86_64:$NGSchoolDir/src/
export PATH=$PATH:$NGSchoolDir/src/augustus/src::$NGSchoolDir/src/augustus-3.0.2/bin

export PATH=$PATH:$NGSchoolDir/src/snpEff

export PATH=$PATH:$NGSchoolDir/src/Bismark:$NGSchoolDir/src/preseq_v2.0:$NGSchoolDir/src/qualimap_v2.2:$NGSchoolDir/src/trim_galore_zip
export PATH=$PATH::$NGSchoolDir/src/convector

export SSPACEBIN=$NGSchoolDir/src/SSPACE/SSPACE_Standard_v3.0.pl

echo ".bashrc loaded!"
