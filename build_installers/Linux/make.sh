#!/bin/sh

### This file will not install myPhyloDB.
### Please download the installer at: http://www.ars.usda.gov/services/software/download.htm?softwareid=472

cd $HOME/PycharmProjects/myPhyloDB/
pyinstaller -D $HOME/PycharmProjects/myPhyloDB/serve-linux.spec
cd $HOME/PycharmProjects/myPhyloDB/dist/
tar -zcvf myPhyloDB.tar.gz myPhyloDB/*
mv myPhyloDB.tar.gz $HOME/PycharmProjects/myPhyloDB/build_installers/Linux/myPhyloDB.tar.gz
cd $HOME/PycharmProjects/myPhyloDB/build_installers/Linux
chmod +x ./install.sh
$HOME/megastep-makeself/makeself.sh . myPhyloDB_v.1.1.2_Linux_x64_install.sh "myPhyloDB installer..." ./install.sh
