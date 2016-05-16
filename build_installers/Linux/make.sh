#!/bin/sh
### This file is used to make a new myPhyloDB installer
#    it will not install myPhyloDB.

### If you are looking for a pre-built installer
# Please download the installer at: http://www.ars.usda.gov/services/software/download.htm?softwareid=472

### To run this file type the following in your terminal
# cd $HOME/PycharmProjects/myPhyloDB/build_installers/Linux
# workon myphylodb
# sh make.sh

cd $HOME/PycharmProjects/myPhyloDB
pyinstaller -D $HOME/PycharmProjects/myPhyloDB/serve-linux.spec
deactivate
cd $HOME/PycharmProjects/myPhyloDB/dist/
tar -zcvf myPhyloDB.tar.gz myPhyloDB/*
mv myPhyloDB.tar.gz $HOME/PycharmProjects/myPhyloDB/build_installers/Linux/myPhyloDB.tar.gz
cd $HOME/PycharmProjects/myPhyloDB/build_installers/Linux
chmod +x ./install.sh
$HOME/megastep-makeself/makeself.sh . myPhyloDB_v.1.2.0_Linux_x64_install.sh "myPhyloDB installer..." ./install.sh
