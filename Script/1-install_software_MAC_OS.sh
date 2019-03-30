# download and install QIIME2 following to web site;
wget https://data.qiime2.org/distro/core/qiime2-2019.1-py36-osx-conda.yml
conda env create -n qiime2-2019.1 --file qiime2-2019.1-py36-osx-conda.yml
# OPTIONAL CLEANUP
rm qiime2-2019.1-py36-osx-conda.yml
