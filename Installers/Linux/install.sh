#! /usr/bin/sh

tar -zxf myPhyloDB.tar.gz -C $HOME
mkdir -p $HOME/.icons
cp $HOME/myPhyloDB/media/images/database_2_48.png $HOME/.icons/database_2_48.png
cp myPhyloDB.desktop $HOME/myPhyloDB/myPhyloDB.desktop
chmod +x $HOME/myPhyloDB/myPhyloDB.desktop

echo Done!
