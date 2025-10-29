#!/bin/bash

chmod +x dScaff.sh
chmod +x contigs_indexing.py

cd ~/Downloads
sudo apt install -y ncbi-blast+

sudo apt install -y git
sudo apt install -y seqtk

#git clone https://github.com/lh3/seqtk.git
#cd seqtk; make

cd ~/Downloads
sudo apt install -y r-base

if ! command -v python3 &> /dev/null; then
    echo "Python3 is not installed. Installing..."
    sudo apt update
    sudo apt install python3 -y
else
    echo "Python3 is already installed: $(python3 --version)"
fi



