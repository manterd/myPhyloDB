#! /usr/bin/sh

echo "Checking if a current version of myPhyloDB exists..."

if [-d $HOME/myPhyloDB]
    then
        echo -n "Would you like to keep your current database files [y/n]?"
        read response
    else
        response=y
fi

echo response

cp $HOME/myPhyloDB/dbMicrobe $HOME/myPhyloDB_temp/dbMicrobe
cp -r $HOME/myPhyloDB/uploads $HOME/myPhyloDB_temp/uploads
tar -zxf myPhyloDB.tar.gz -C $HOME
mv ./dbMicrobe $HOME/dbMicrobe
mv ./uploads $HOME/uploads

mkdir -p $HOME/.icons
cp $HOME/myPhyloDB/media/images/database_2_48.png $HOME/.icons/database_2_48.png
cp myPhyloDB.desktop $HOME/Desktop/myPhyloDB.desktop
chmod +x $HOME/Desktop/myPhyloDB.desktop

echo Done!
