#! /usr/bin/sh

mv $HOME/dbMicrobe .
mv $HOME/uploads .
tar -zxf myPhyloDB.tar.gz -C $HOME
mkdir -p $HOME/.icons
cp $HOME/myPhyloDB/media/images/database_2_48.png $HOME/.icons/database_2_48.png
cp myPhyloDB.desktop $HOME/Desktop/myPhyloDB.desktop
chmod +x $HOME/Desktop/myPhyloDB.desktop

echo Done!
