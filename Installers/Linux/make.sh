#!/bin/sh

cd $HOME/PycharmProjects/myPhyloDB/
pyinstaller -D $HOME/PycharmProjects/myPhyloDB/serve-linux.spec
cd $HOME/PycharmProjects/myPhyloDB/dist/
tar -zcvf myPhyloDB.tar.gz myPhyloDB/*
mv myPhyloDB.tar.gz $HOME/PycharmProjects/myPhyloDB/Installers/Linux/myPhyloDB.tar.gz
cd $HOME/PycharmProjects/myPhyloDB/Installers/Linux
chmod +x ./install.sh
$HOME/megastep-makeself-be1c982/makeself.sh . myPhyloDB_1.0_Linux_x64_install.sh "myPhyloDB installer..." ./install.sh
