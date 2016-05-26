#!/usr/bin/env bash

# Aim: download and install packages in my home
# Copyright (C) 2014-2016 Institut National de la Recherche Agronomique
# License: GPL-3+
# Author: Timothée Flutre
# Versioning: https://github.com/timflutre/quantgen

set -e

# list of programs in alphabetical order
declare -a progs=("art" "artfastqgen" "autoconf" "automake" "bedtools" "biobambam" "blup_gen_snp" "bsfg" "bwa" "carthagene" "cutadapt" "deindexer" "dmu" "dnemulator" "dwgsim" "eagle" "ea-utils" "eigensoft" "emacs" "epcr" "eqtlbma" "ess" "fastqc" "forqs" "gbs-barcode-splitter" "gemma" "gs3" "gsl" "htslib" "help2man" "igv" "inphap" "insilicut" "latex2html" "ldso" "libtool" "lsof" "mapmaker" "ms" "mstrat" "openbugs" "patman" "platypus" "primer3" "polymode" "R" "rar" "repet" "rpy2" "quantinemo" "samtools" "scilab" "sickle" "smart" "southgreen_utils" "stacks" "tabula" "tar" "tedna" "texinfo" "texlive" "tm" "tmap" "trim-galore" "trimmomatic" "ubd" "wgsim" "xclip" "zlib")

if [ "$#" -ne 1 ]; then
    echo "ERROR: need to provide a program name as parameter"
    exit 1
fi

# http://stackoverflow.com/a/8574392/597069
containsElement () {
    local e
    for e in "${@:2}"; do [[ "$e" == "$1" ]] && return 0; done
    return 1
}

if ! containsElement "$1" "${progs[@]}"; then
    echo -e "WARNING: '$1' is unknown to me"
    exit 0
fi

date

if [ "$1" == "art" ]; then
    mkdir -p $1
    cd $1
    wget http://www.niehs.nih.gov/research/resources/assets/docs/artalllinux64bin_bptargz.gz
    tar -xzvf artalllinux64bin_bptargz.gz
    cd Linux64
    cp 454_readprofile_art aln2bed.pl art_454 art_illumina art_SOLiD $HOME/bin/
fi

if [ "$1" == "artfastqgen" ]; then
    mkdir -p $1
    cd $1
    wget http://sourceforge.net/projects/artfastqgen/files/latest/download
    unzip ArtificialFastqGenerator.zip
fi

if [ "$1" == "autoconf" ]; then
    mkdir -p $1
    cd $1
    wget ftp://ftp.gnu.org/gnu/autoconf/autoconf-2.69.tar.gz
    tar -xzvf autoconf-2.69.tar.gz
    cd autoconf-2.69
    ./configure --prefix=$HOME
    make
    # make check
    make install
fi

if [ "$1" == "automake" ]; then
    mkdir -p $1
    cd $1
    wget ftp://ftp.gnu.org/gnu/automake/automake-1.13.1.tar.gz
    tar -xzvf automake-1.13.1.tar.gz
    cd automake-1.13.1
    ./configure --prefix=$HOME
    make
    # make check
    make install
fi

if [ "$1" == "bedtools" ]; then
    mkdir -p $1
    cd $1
    wget --no-check-certificate -O bedtools-2.19.1.tar.gz https://github.com/arq5x/bedtools2/releases/download/v2.19.1/bedtools-2.19.1.tar.gz
    tar -xzvf bedtools-2.19.1.tar.gz
    cd bedtools2-2.19.1
    make
    cp bin/* $HOME/bin/
fi

if [ "$1" == "biobambam" ]; then
    mkdir -p $1
    cd $1
    wget https://github.com/gt1/biobambam2/releases/download/2.0.9-release-20150619154907/biobambam2-2.0.9-release-20150619154907-x86_64-etch-linux-gnu.tar.gz
    tar -xzvf biobambam2-2.0.9-release-20150619154907-x86_64-etch-linux-gnu.tar.gz
    cd biobambam2-2.0.9-release-20150619154907-x86_64-etch-linux-gnu/
    cp bin/bamsormadup \
       bin/bamcollate2 \
       bin/bammarkduplicates2 \
       bin/bamtofastq \
       $HOME/bin/
fi

if [ "$1" == "blup_gen_snp" ]; then
    mkdir -p $1
    cd $1
    wget http://snp.toulouse.inra.fr/~alegarra/progs_genom_sel.tar.gz
    tar -xzvf progs_genom_sel.tar.gz
    cd progs_genom_sel
    cp blup_gen blup_snp $HOME/bin
fi

if [ "$1" == "bsfg" ]; then
    mkdir -p $1
    cd $1
    wget http://www.stat.duke.edu/~sayan/bfgr/BSF-G.zip
    unzip BSF-G.zip
fi

if [ "$1" == "bwa" ]; then
    mkdir -p $1
    cd $1
    wget http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.12.tar.bz2/download
    tar -xvf bwa-0.7.12.tar.bz2
    cd bwa-0.7.12/
    make
    cp bwa $HOME/bin
fi

if [ "$1" == "carthagene" ]; then
    mkdir -p $1
    cd $1
    # git clone http://mulcyber.toulouse.inra.fr/anonscm/git/carthagene/carthagene.git
    git clone -b v1.3 http://mulcyber.toulouse.inra.fr/anonscm/git/carthagene/carthagene.git
    cd carthagene
    CGROOT=$(pwd)
    wget http://www7.inra.fr/mia/T/CarthaGene/Download/CarthaGene-server.r
    cd ..; mkdir build; cd build
    cmake $CGROOT -DWITH_LKH=ON -DWITH_FRAMEWORK=ON -DCMAKE_INSTALL_PREFIX=$HOME -DCMAKE_BUILD_TYPE=Release #-DGENERATE_DOC=ON
    sed -i 's:/\*#include "Marker.h"\*/:#include "Marker.h":' ../carthagene/calcul/CartaGene.h
    make
    cp bin/carthagene $HOME/bin
    cp -r share/* $HOME/share
    # get the doc from the official release
    cd ..
    wget https://mulcyber.toulouse.inra.fr/frs/download.php/1196/carthagene-nonfree-1.3.beta-Linux-x86_64.sh
    sh ./carthagene-nonfree-1.3.beta-Linux-x86_64.sh <<EOF
y
N
EOF
fi

if [ "$1" == "cutadapt" ]; then
    mkdir -p $1
    cd $1
    # wget --no-check-certificate -O master.zip https://github.com/marcelm/cutadapt/archive/master.zip
    # unzip master.zip
    # cd cutadapt-master/
    wget --no-check-certificate -O cutadapt-1.8.tar.gz https://github.com/marcelm/cutadapt/archive/v1.8.tar.gz
    tar -xvf cutadapt-1.8.tar.gz
    cd cutadapt-1.8/
    python setup.py install --prefix=$HOME
    echo "be sure to include ${HOME}/lib/... in your PYTHONPATH"
fi

if [ "$1" == "deindexer" ]; then
    mkdir -p $1
    cd $1
    wget --no-check-certificate -O master.zip https://github.com/ws6/deindexer/archive/master.zip
    unzip master.zip
    cd deindexer-master
    make
    cp deindexer $HOME/bin
fi

if [ "$1" == "dmu" ]; then
    mkdir -p $1
    cd $1
    wget http://dmu.agrsci.dk/DMU/Linux/Current/dmuv6-R5.2-EM64T-build-2014-10-23.tar.gz
    tar -xzvf dmuv6-R5.2-EM64T-build-2014-10-23.tar.gz
    cd dmuv6/R5.2-EM64T/; mkdir -p doc; cd doc; \
        wget http://dmu.agrsci.dk/DMU/Doc/Current/dmuv6_guide.5.2.pdf
    cd ..; cp bin/* $HOME/bin
fi

if [ "$1" == "dnemulator" ]; then
    mkdir -p $1
    cd $1
    wget http://www.cbrc.jp/dnemulator/dnemulator-16.zip
    unzip dnemulator-16.zip
    make install prefix=$HOME
fi

if [ "$1" == "dwgsim" ]; then
    mkdir -p $1
    cd $1
    # wget http://sourceforge.net/projects/dnaa/files/latest/download?source=files
    # tar -xzvf dwgsim-0.1.10.tar.gz
    # cd dwgsim-0.1.10
    # wget http://sourceforge.net/projects/samtools/files/samtools/0.1.7/samtools-0.1.7a.tar.bz2/download
    # tar -xjvf samtools-0.1.7a.tar.bz2
    # ln -s samtools-0.1.7a samtools
    # wget http://sourceforge.net/projects/samtools/files/samtools/0.1.8/samtools-0.1.8.tar.bz2/download
    # tar -xjvf samtools-0.1.8.tar.bz2
    # ln -s samtools-0.1.8 samtools
    # wget http://sourceforge.net/projects/samtools/files/latest/download?source=files
    # tar -xjvf samtools-0.1.19.tar.bz2
    # ln -s samtools-0.1.19 samtools
    # make
    git clone git@github.com:nh13/DWGSIM.git
    cd DWGSIM
    git submodule init
    git submodule update
    make
    cp dwgsim dwgsim_eval dwgims_pileup_eval.pl $HOME/bin
fi

if [ "$1" == "eagle" ]; then
    mkdir -p $1
    cd $1
    wget --no-check-certificate -O master.zip https://github.com/sequencing/EAGLE/archive/master.zip
    unzip master.zip
    cd EAGLE-master/
    ./src/configure --prefix=$HOME # requires boost
fi

if [ "$1" == "ea-utils" ]; then
    wget 
    tar -xzvf ea-utils.1.1.2-806.tar.gz
    cd ea-utils.1.1.2-806
    PREFIX=$HOME make install
fi

if [ "$1" == "eigensoft" ]; then
    wget ftp://pricelab:pricelab@ftp.broadinstitute.org/EIGENSOFT/EIG6.0beta.tar.gz
    tar -xzvf EIG6.0beta.tar.gz
    cd EIG6.0beta
    cp bin/* $HOME/bin
fi

if [ "$1" == "emacs" ]; then
    mkdir -p $1
    cd $1
    wget --timestamping http://gnu.mirrors.hoobly.com/gnu/emacs/emacs-24.5.tar.gz
    tar -xzvf emacs-24.5.tar.gz
    cd emacs-24.5
    ./configure --prefix=$HOME #--with-x-toolkit=no --with-xpm=no --with-jpeg=no --with-gif=no --with-tiff=no
    make
    make install
fi

if [ "$1" == "epcr" ]; then
    mkdir -p $1
    cd $1
    wget ftp://ftp.ncbi.nlm.nih.gov/pub/schuler/e-PCR/e-PCR-2.3.12-1-src.tar.gz
    tar -xzvf e-PCR-2.3.12-1-src.tar.gz
    cd e-PCR-2.3.12
    
fi

if [ "$1" == "eqtlbma" ]; then
    mkdir -p $1
    cd $1
    # wget 
    # tar -xzvf 
    git clone https://github.com/timflutre/eqtlbma.git
    cd eqtlbma
    autoreconf --install --force --symlink
    ./configure --prefix=$HOME
    make
    # make check
    make install
fi

if [ "$1" == "ess" ]; then
    mkdir -p $1
    cd $1
    wget -O ess-15.09-2.tgz http://ess.r-project.org/downloads/ess/ess-15.09-2.tgz
    tar -xzvf ess-15.09-2.tgz
    cd ess-15.09-2
    make
    echo "read the manual to complete the install"
    echo "http://ess.r-project.org/Manual/ess.html#Installation"
fi

if [ "$1" == "fastqc" ]; then
    mkdir -p $1
    cd $1
    wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.2.zip
    unzip fastqc_v0.11.2.zip
    cd FastQC
    chmod +x fastqc
    echo "add "$(pwd)" to your PATH"
fi

if [ "$1" == "forqs" ]; then
    mkdir -p $1
    cd $1
    wget --no-check-certificate https://bitbucket.org/dkessner/forqs/downloads/forqs_linux_131213_104438.zip
    unzip forqs_linux_131213_104438.zip
    cd forqs_linux_131213_104438
    cp bin/* $HOME/bin
fi

if [ "$1" == "gbs-barcode-splitter" ]; then
    mkdir -p $1
    cd $1
    wget http://sourceforge.net/projects/gbsbarcode/files/GBS_barcode.pl/download
    wget http://sourceforge.net/projects/gbsbarcode/files/apeki.txt/download
    wget http://sourceforge.net/projects/gbsbarcode/files/GBS_barcode_template.txt/download
    wget http://sourceforge.net/projects/gbsbarcode/files/readme.txt/download
fi

if [ "$1" == "gemma" ]; then
    mkdir -p $1
    cd $1
    wget http://home.uchicago.edu/xz7/software/gemma-0.94.tar.gz
    mkdir gemma-0.94; cd gemma-0.94; ln -s ../gemma-0.94.tar.gz .
    tar -xzvf gemma-0.94.tar.gz
    cd 
    make
    cp bin/gemma $HOME/bin
    # git clone git@github.com:xiangzhou/GEMMA.git
    # put FORCE_DYNAMIC = 1
    # cp bin/gemma $HOME/bin
fi

if [ "$1" == "gs3" ]; then
    mkdir -p $1
    cd $1
    wget http://snp.toulouse.inra.fr/~alegarra/gs3dist.tar.gz
    tar -xzvf gs3dist.tar.gz
    cd gs3dist
    autogen.sh
    # ./configure --prefix=$HOME
    # make
    # make install
    wget http://snp.toulouse.inra.fr/~alegarra/gs3_linux64bit_executable
    chmod +x gs3_linux64bit_executable
    cp gs3_linux64bit_executable $HOME/bin
fi

if [ "$1" == "gsl" ]; then
    mkdir -p $1
    cd $1
    wget ftp://ftp.gnu.org/gnu/gsl/gsl-1.16.tar.gz
    tar -xzvf gsl-1.16.tar.gz
    cd gsl-1.16
    ./configure --prefix=$HOME
    make
    # make check
    make install
fi

if [ "$1" == "htslib" ]; then
    mkdir -p $1
    cd $1
    wget -O htslib-1.2.1.tar.bz2 --no-check-certificate https://github.com/samtools/htslib/releases/download/1.2.1/htslib-1.2.1.tar.bz2
    tar -xvf htslib-1.2.1.tar.bz2
    cd htslib-1.2.1/
    ./configure --prefix=$HOME
    make
    make install
fi

if [ "$1" == "help2man" ]; then
    mkdir -p $1
    cd $1
    wget http://mirror.ibcp.fr/pub/gnu/help2man/help2man-1.43.3.tar.gz
    tar -xzvf help2man-1.43.3.tar.gz
    cd help2man-1.43.3
    ./configure --prefix=$HOME
    make
    make install
fi

if [ "$1" == "igv" ]; then
    mkdir -p $1
    cd $1
    # wget https://github.com/igvteam/igv/archive/v2.3.57.zip
    # unzip v2.3.57.zip
    # cd igv-2.3.57/
    wget http://data.broadinstitute.org/igv/projects/downloads/IGV_2.3.57.zip
    unzip IGV_2.3.57.zip
    cd IGV_2.3.57/
    cp igv.jar igv.sh $HOME
fi

if [ "$1" == "inphap" ]; then
    mkdir -p $1
    cd $1
    wget http://it.informatik.uni-tuebingen.de/software/inPhap/inPhap.jar
    wget http://it.informatik.uni-tuebingen.de/software/inPhap/inPHAP_Example_Files.zip
fi

if [ "$1" == "insilicut" ]; then
    mkdir -p $1
    cd $1
    wget -O insilicut-1.1.2.tar.gz https://github.com/timflutre/insilicut/archive/v1.1.2.tar.gz
    tar -xzvf insilicut-1.1.2.tar.gz
    cd insilicut-1.1.2/
    make
    # make check
    make install
fi

if [ "$1" == "latex2html" ]; then
    mkdir -p $1
    cd $1
    wget http://mirrors.ctan.org/support/latex2html/latex2html-2012.tgz
    tar -xzvf latex2html-2012.tgz
    cd latex2html-2012
    ./configure --prefix=$HOME
    make
    make install
fi

if [ "$1" == "ldso" ]; then
    mkdir -p $1
    cd $1
    wget https://qgsp.jouy.inra.fr/archives/LDSO/LDSO_v1.02.rar
    mkdir LDSO_v1.02; cd LDSO_v1.02; ln -s ../LDSO_v1.02.rar .
    rar e LDSO_v1.02.rar <<EOF
A
EOF
    chmod +x Ldso_v1.02.f90
    cp Ldso_v1.02.f90 $HOME/bin
fi

if [ "$1" == "libtool" ]; then
    mkdir -p $1
    cd $1
    wget ftp://ftp.gnu.org/gnu/libtool/libtool-2.4.2.tar.gz
    tar -xzvf libtool-2.4.2.tar.gz
    cd libtool-2.4.2
    ./configure --prefix=$HOME
    make
    # make check
    make install
fi

if [ "$1" == "lsof" ]; then
    mkdir -p $1
    cd $1
    wget ftp://lsof.itap.purdue.edu/pub/tools/unix/lsof/lsof_4.87.tar.gz
    tar -xzvf lsof_4.87.tar.gz
    cd lsof_4.87
    tar xf lsof_4.87_src.tar
    cd lsof_4.87_src
    ./Configure linux
    make
fi

if [ "$1" == "mapmaker" ]; then
    mkdir -p $1
    cd $1
    wget http://www.broadinstitute.org/ftp/distribution/software/mapmaker3/mapm3-source.tar.Z
    tar -xzvf mapm3-source.tar.Z
    sed -i 's/SYS= -D_SYS_WATCOM/SYS= -D_SYS_WATCOM/' Makefile
    make
    echo "see http://www.unix.com/programming/169884-compiling-old-c-program-linux.html"
    echo "see also http://www.latitudecartography.co.uk/forum/topic.asp?TOPIC_ID=2033"
fi

if [ "$1" == "ms" ]; then
    mkdir -p $1
    cd $1
    wget https://webshare.uchicago.edu/users/rhudson1/Public/ms.folder/ms.tar.gz
    tar -xzvf ms.tar.gz
    cd msdir/
    gcc -O3 -o ms ms.c streec.c rand1.c -lm # uses drand48()
    cp ms $HOME/bin
fi

if [ "$1" == "mstrat" ]; then
    mkdir -p $1
    cd $1
    wget http://www.ensam.inra.fr/gap/MSTRAT/download/MStrat-v4-unix.zip
    unzip MStrat-v4-unix.zip
    cd MStrat-v4-unix
    gcc -o ether Etherv4.c -lm
fi

if [ "$1" == "openbugs" ]; then
    mkdir -p $1
    cd $1
    wget http://www.openbugs.net/w/OpenBUGS_3_2_3?action=AttachFile&do=get&target=OpenBUGS-3.2.3.tar.gz
    tar -xzvf OpenBUGS-3.2.3.tar.gz
    cd OpenBUGS-3.2.3/
    ./configure --prefix=$HOME
    make
    make install
fi

if [ "$1" == "patman" ]; then
    mkdir -p $1
    cd $1
    wget --no-check-certificate https://bioinf.eva.mpg.de/patman/patman-1.2.2.tar.gz
    tar -xzvf patman-1.2.2.tar.gz
    cd patman-1.2.2
    make
    make install
fi

if [ "$1" == "platypus" ]; then
    mkdir -p $1
    cd $1
    # need to register: http://www.well.ox.ac.uk/platypus
    wget http://www.well.ox.ac.uk/bioinformatics/Software/Platypus-latest.tgz
    tar -xvf Platypus-latest.tgz
    cd Platypus_0.8.1/
    ./buildPlatypus.sh
    python setup.py install --home=$HOME
    echo -e '#!/usr/bin/env python\n' > tmp
    cat Platypus.py >> tmp
    mv tmp Platypus.py
    cp Platypus.py $HOME/bin
fi

if [ "$1" == "primer3" ]; then
    mkdir -p $1
    cd $1
    wget http://downloads.sourceforge.net/project/primer3/primer3/2.3.6/primer3-src-2.3.6.tar.gz
    tar -xzvf primer3-src-2.3.6.tar.gz
    cd tar -xzvf primer3-2.3.6/
    cd src
    make
    cp primer3_core $HOME/bin
fi

if [ "$1" == "polymode" ]; then
    mkdir -p $1
    cd $1
    git clone https://github.com/vitoshka/polymode.git
fi

if [ "$1" == "quantinemo" ]; then
    mkdir -p $1
    cd $1
    wget http://www2.unil.ch/popgen/softwares/quantinemo/quantinemo_files/quantinemo_src.zip
    unzip quantinemo_src.zip
    cd quantinemo_src
    make
    cp quantiNemo $HOME/bin
fi

if [ "$1" == "samtools" ]; then
    mkdir -p $1
    cd $1
    # wget http://sourceforge.net/projects/samtools/files/samtools/1.0/samtools-1.0.tar.bz2/download
    tar -xjvf samtools-1.0.tar.bz2
    wget http://sourceforge.net/projects/samtools/files/samtools/1.1/samtools-1.1.tar.bz2/download
    tar -xjvf samtools-1.1.tar.bz2
    cd samtools-1.1/
    make
    make prefix=$HOME install
    cd htslib-1.1/; make bgzip; make tabix; cp bgzip tabix $HOME/bin
fi

if [ "$1" == "scilab" ]; then
    mkdir -p $1
    cd $1
    wget http://www.scilab.org/download/5.5.1/scilab-5.5.1.bin.linux-x86_64.tar.gz
    tar -xzvf scilab-5.5.1.bin.linux-x86_64.tar.gz
fi

if [ "$1" == "sickle" ]; then
    mkdir -p $1
    cd $1
    wget --no-check-certificate -O sickle-master.zip https://github.com/najoshi/sickle/archive/master.zip
    unzip sickle-master.zip
    cd sickle-master
    ./configure --prefix=$HOME
    make
    cp sickle $HOME/bin
fi

if [ "$1" == "smart" ]; then
    mkdir -p $1
    cd $1
    wget https://urgi.versailles.inra.fr/download/s-mart/s-mart-1.1.4.zip
    unzip s-mart-1.1.4.zip
    # cd 
    # cp sickle $HOME/bin
fi

if [ "$1" == "southgreen_utils" ]; then
    mkdir -p $1
    cd $1
    wget --no-check-certificate -O master.zip https://github.com/SouthGreenPlatform/utils/archive/master.zip
    unzip master.zip
    cd utils-master
    chmod +x transpose_annotation/transpose_annotation.pl
    cp transpose_annotation/transpose_annotation.pl $HOME/bin
fi

if [ "$1" == "stacks" ]; then
    mkdir -p $1
    cd $1
    wget http://creskolab.uoregon.edu/stacks/source/stacks-1.27.tar.gz
    tar -xzvf stacks-1.27.tar.gz
    cd stacks-1.27
    ./configure --prefix=$HOME --enable-sparsehash --enable-bam
    make
    make install
fi

if [ "$1" == "R" ]; then
    mkdir -p $1
    cd $1
    wget http://cran.rstudio.com/src/base/R-3/R-3.3.0.tar.gz
    tar -xvf R-3.3.0.tar.gz
    cd R-3.3.0
    ./configure --prefix=$HOME
    make
    make install
fi

if [ "$1" == "rar" ]; then
    mkdir -p $1
    cd $1
    wget http://www.rarlab.com/rar/rarlinux-x64-5.0.1.tar.gz
    tar -xzvf rarlinux-x64-5.0.1.tar.gz
    cd rar
    # make install PREFIX=$HOME
    cp rar_static $HOME/bin/rar
fi

if [ "$1" == "repet" ]; then
    mkdir -p $1
    cd $1
    wget https://urgi.versailles.inra.fr/download/repet/REPET_linux-x64-2.2.tar.gz
    tar -xzvf REPET_linux-x64-2.2.tar.gz
    cd REPET_linux-x64-2.2
    cp bin/* $HOME/bin
fi

if [ "$1" == "rpy2" ]; then
    mkdir -p $1
    cd $1
    wget https://pypi.python.org/packages/source/r/rpy2/rpy2-2.6.0.tar.gz#md5=679898fbc832d4f05a5efcf1a7eb1a68
    tar -xzvf rpy2-2.6.0.tar.gz
    cd rpy2-2.6.0
    python setup.py install --user
fi

if [ "$1" == "tabula" ]; then
    mkdir -p $1
    cd $1
    wget --no-check-certificate -O tabula-jar-0.9.6.zip https://github.com/tabulapdf/tabula/releases/download/v0.9.6/tabula-jar-0.9.6.zip
    unzip tabula-jar-0.9.6.zip
    cd tabula
    echo '#!/usr/bin/env bash\njava -Dfile.encoding=utf-8 -Xms256M -Xmx1024M -jar tabula.jar' > tabula
    chmod +x tabula
    cp tabula tabula.jar $HOME/bin
fi

if [ "$1" == "tar" ]; then
    mkdir -p $1
    cd $1
    wget ftp://ftp.gnu.org/gnu/tar/tar-1.27.tar.gz
    tar -xzvf tar-1.27.tar.gz
    cd tar-1.27
    ./configure --prefix=$HOME
    make
    # make check
    make install
fi

if [ "$1" == "texlive" ]; then
    mkdir -p $1
    cd $1
    wget http://mirror.ctan.org/systems/texlive/tlnet/install-tl-unx.tar.gz
    tar -xzvf install-tl-unx.tar.gz
    cd install-tl-20140918
    ./install-tl
fi

if [ "$1" == "tedna" ]; then
    mkdir -p $1
    cd $1
    wget https://urgi.versailles.inra.fr/content/download/3481/29402/file/tedna_1.2.2.tar.gz
    tar -xzvf tedna_1.2.2.tar.gz
    make
    cp tedna $HOME/bin
    # cp Evaluator/Evaluator.py $HOME/bin # requires S-MART
fi

if [ "$1" == "texinfo" ]; then
    mkdir -p $1
    cd $1
    wget http://ftp.gnu.org/gnu/texinfo/texinfo-5.2.tar.gz
    tar -xzvf texinfo-5.2.tar.gz
    cd texinfo-5.2
    ./configure --prefix=$HOME
    make
    # make check
    make install
fi

if [ "$1" == "tm" ]; then
    mkdir -p $1
    cd $1
    wget http://snp.toulouse.inra.fr/~alegarra/TMdist.tar.gz
    tar -xzvf TMdist.tar.gz
    cd TMdist
    wget http://snp.toulouse.inra.fr/~alegarra/manualtm.pdf
    cp tm $HOME/bin
fi

if [ "$1" == "tmap" ]; then
    mkdir -p $1
    cd $1
    wget http://users.math.yale.edu/~dc597/tmap/tmap-1.1.tar.gz
    tar -xzvf tmap-1.1.tar.gz
    cd tmap
    ./configure --prefix=$HOME
    make
    make install
fi

if [ "$1" == "trim-galore" ]; then
    mkdir -p $1
    cd $1
    wget http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/trim_galore_v0.3.7.zip
    unzip trim_galore_v0.3.7.zip
    cd trim_galore_zip/
    cp trim_galore $HOME/bin
fi

if [ "$1" == "trimmomatic" ]; then
    mkdir -p $1
    cd $1
    # wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-Src-0.32.zip
    # unzip Trimmomatic-Src-0.32.zip
    wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.33.zip
    unzip Trimmomatic-0.33.zip
    cd Trimmomatic-0.33/
    cp trimmomatic-0.33.jar $HOME/bin
fi

if [ "$1" == "ubd" ]; then
    mkdir -p $1
    cd $1
    wget --no-check-certificate -O master.zip https://github.com/pelinakan/UBD/archive/master.zip
    unzip master.zip
    cd UBD-master
    # make
    cp bin/linux/* $HOME/bin
fi

if [ "$1" == "wgsim" ]; then
    mkdir -p $1
    cd $1
    wget --no-check-certificate -O master.zip https://github.com/lh3/wgsim/archive/master.zip
    unzip master.zip
    cd wgsim-master/
    gcc -g -O2 -Wall -o wgsim wgsim.c -lz -lm
    cp wgsim wgsim_eval.pl $HOME
fi

if [ "$1" == "xclip" ]; then
    mkdir -p $1
    cd $1
    wget --no-check-certificate http://sourceforge.net/projects/xclip/files/latest/download
    tar -xzvf xclip-0.12.tar.gz
    cd xclip-0.12
    ./configure --prefix=$HOME
    make
    make install
fi

if [ "$1" == "zlib" ]; then
    mkdir -p $1
    cd $1
    wget http://zlib.net/zlib-1.2.8.tar.gz
    tar -xzvf zlib-1.2.8.tar.gz
    cd zlib-1.2.8
    ./configure --prefix=$HOME
    make
    # make check
    make install
fi

date
